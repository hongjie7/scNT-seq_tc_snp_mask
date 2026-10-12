#!/usr/bin/env python3
"""
plot_tc_reads.py

scNT-seq (Qiu et al. 2020, Fig. 1c)-style plot: every unique transcript (UMI)
of one gene in one cell, with each covered uridine shown as a grey dot and each
uridine carrying a T-to-C substitution shown as an "X" coloured by how many
reads of that UMI support the substitution.

Optionally overlays a SNP mask (e.g. your WGS T>C / A>G CSV) so you can see
whether genetic variants are masked correctly.

IMPORTANT: use the ORIGINAL (pre-consensus) BAM, e.g. donor_bams/<donor>.bam or
STARsolo Aligned.sortedByCoord.out.bam. In the dynast consensus BAM each UMI is
collapsed to one read, so the "coverage" colour would always be 1.

Strand logic:
    gene on + strand: uridine = genomic T, conversion = T>C
    gene on - strand: uridine = genomic A, conversion = A>G

Examples:
    # list the cells with most UMIs for a gene, then plot the top one
    python plot_tc_reads.py --bam donor_bams/DONOR.bam \
        --fasta GRCh38.p14.genome.fa --gtf gencode.v47.annotation.gtf \
        --gene FOS --out FOS_topcell.pdf

    # a specific cell, with the SNP mask overlaid
    python plot_tc_reads.py --bam donor_bams/DONOR.bam \
        --fasta GRCh38.p14.genome.fa --gtf gencode.v47.annotation.gtf \
        --gene FOS --cell AAACCTGAGAAACCAT \
        --snp-csv WGS_SNP_mask_0-based/DONOR_snps.csv --snp-coord 0based \
        --out FOS_AAACCTGAGAAACCAT.pdf
"""
import argparse
import gzip
import sys
from collections import defaultdict, Counter

import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D


# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------
def find_gene(gtf, gene, bam_contigs):
    """Return (chrom, start0, end, strand, gene_id) for gene_name == gene."""
    opener = gzip.open if gtf.endswith(".gz") else open
    key = f'gene_name "{gene}";'
    hits = []
    with opener(gtf, "rt") as f:
        for line in f:
            if line.startswith("#") or key not in line:
                continue
            p = line.rstrip("\n").split("\t")
            if p[2] != "gene":
                continue
            gid = p[8].split('gene_id "')[1].split('"')[0]
            hits.append((p[0], int(p[3]) - 1, int(p[4]), p[6], gid))
    if not hits:
        sys.exit(f"ERROR: gene_name '{gene}' not found in GTF")
    # prefer contigs present in the BAM and primary chromosomes (no '_' patches)
    hits.sort(key=lambda h: (h[0] not in bam_contigs, "_" in h[0]))
    if hits[0][0] not in bam_contigs:
        sys.exit(f"ERROR: none of {[h[0] for h in hits]} are contigs in the BAM")
    return hits[0]


def load_snps(snp_csv, chrom, coord):
    """Load masked positions on `chrom` as a set of 0-based positions.
    Works with or without a 'contig,genome_i,...' header."""
    pos = set()
    offset = 1 if coord == "1based" else 0
    with open(snp_csv) as f:
        for i, line in enumerate(f):
            p = line.strip().split(",")
            if len(p) < 2:
                continue
            if i == 0 and not p[1].isdigit():   # header line
                continue
            if p[0] == chrom:
                pos.add(int(p[1]) - offset)
    return pos


def norm_bc(bc):
    return bc[:-2] if bc.endswith("-1") else bc


def strand_ok(read, gene_strand, lib_strand):
    if lib_strand == "unstranded":
        return True
    sense = read.is_reverse == (gene_strand == "-")
    return sense if lib_strand == "forward" else not sense


def iter_gene_reads(bam, chrom, start, end, gene_strand, lib_strand,
                    cb_tag, ub_tag, mapq):
    for r in bam.fetch(chrom, start, end):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if r.mapping_quality < mapq:
            continue
        if not (r.has_tag(cb_tag) and r.has_tag(ub_tag)):
            continue
        if not strand_ok(r, gene_strand, lib_strand):
            continue
        yield r


# ----------------------------------------------------------------------------
# Core: per-UMI uridine coverage and T-to-C counts
# ----------------------------------------------------------------------------
def collect(bam, fasta, chrom, start, end, gene_strand, cell, args):
    ref_u = "T" if gene_strand == "+" else "A"      # uridine on genome
    conv = "C" if gene_strand == "+" else "G"       # converted base on genome
    covered = defaultdict(set)                      # umi -> {pos0}
    conv_n = defaultdict(Counter)                   # umi -> {pos0: n reads}
    n_reads = Counter()

    # dynast 1.0.1 consensus.bam can contain a single-read molecule twice: the
    # original read (no RN tag) and its consensus copy (RN:i:1). If a UMI has any
    # consensus record (RN tag), use only those, so no read is counted twice.
    by_umi = defaultdict(list)
    for r in iter_gene_reads(bam, chrom, start, end, gene_strand, args.strand,
                             args.cb_tag, args.ub_tag, args.min_mapq):
        if norm_bc(r.get_tag(args.cb_tag)) == cell:
            by_umi[r.get_tag(args.ub_tag)].append(r)
    reads = []
    for umi, rs in by_umi.items():
        cons = [r for r in rs if r.has_tag("RN")]
        reads.extend((umi, r) for r in (cons or rs))

    for umi, r in reads:
        n_reads[umi] += r.get_tag("RN") if r.has_tag("RN") else 1
        ref = fasta.fetch(chrom, r.reference_start, r.reference_end).upper()
        seq, qual = r.query_sequence, r.query_qualities
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp < start or rp >= end:
                continue
            if ref[rp - r.reference_start] != ref_u:
                continue
            if qual is not None and qual[qp] < args.min_bq:
                continue
            b = seq[qp]
            if b == "N":
                continue
            covered[umi].add(rp)
            if b == conv:
                conv_n[umi][rp] += 1
    return covered, conv_n, n_reads


def rank_cells(bam, chrom, start, end, gene_strand, args, top=10):
    umis = defaultdict(set)
    for r in iter_gene_reads(bam, chrom, start, end, gene_strand, args.strand,
                             args.cb_tag, args.ub_tag, args.min_mapq):
        umis[norm_bc(r.get_tag(args.cb_tag))].add(r.get_tag(args.ub_tag))
    ranked = sorted(umis.items(), key=lambda kv: -len(kv[1]))
    if not ranked:
        sys.exit("ERROR: no reads with CB/UB found for this gene")
    print("Top cells by UMI count for this gene:")
    for bc, u in ranked[:top]:
        print(f"  {bc}\t{len(u)} UMIs")
    return ranked[0][0]


# ----------------------------------------------------------------------------
# QC: recurrent sites (SNP-like) and mask offset check
# ----------------------------------------------------------------------------
def snp_qc(covered, conv_n, snps, chrom, frac, min_umis):
    cov_umis, conv_umis = Counter(), Counter()
    for umi, sites in covered.items():
        for p in sites:
            cov_umis[p] += 1
            if conv_n[umi][p] > 0:
                conv_umis[p] += 1
    rows = []
    for p, n in cov_umis.items():
        if n >= min_umis and conv_umis[p] / n >= frac:
            rows.append((p, conv_umis[p], n))
    print(f"\nRecurrent T-to-C sites (converted in >= {frac:.0%} of >= "
          f"{min_umis} covering UMIs) -- these look genetic, not 4sU:")
    if not rows:
        print("  none")
    for p, k, n in sorted(rows):
        status = "MASKED" if p in snps else "NOT masked"
        if p not in snps and (p - 1 in snps or p + 1 in snps):
            status += "  <-- mask has a site at +/-1 bp: check --snp-coord / dynast offset!"
        print(f"  {chrom}:{p + 1} (1-based)\t{k}/{n} UMIs\t{status}")
    return rows


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
def plot(covered, conv_n, snps, chrom, gene, gene_strand, cell, args):
    def labeled(umi):
        return any(n >= args.min_conv_reads and p not in snps
                   for p, n in conv_n[umi].items())

    umis = [u for u in covered if covered[u]]
    lab = sorted([u for u in umis if labeled(u)], key=lambda u: min(covered[u]))
    unl = sorted([u for u in umis if not labeled(u)], key=lambda u: min(covered[u]))
    order = (lab + unl)[: args.max_umis]
    n = len(order)
    if n == 0:
        sys.exit("ERROR: this cell has no reads covering uridines of the gene")

    allpos = np.array([p for u in order for p in covered[u]]) + 1
    if args.xlim:
        x0, x1 = args.xlim
    else:
        pad = max(20, int(0.02 * (allpos.max() - allpos.min())))
        x0, x1 = allpos.min() - pad, allpos.max() + pad

    cmap = LinearSegmentedColormap.from_list("bpr", ["#0000ff", "#8000c0", "#ff2020"])
    maxc = max([c for u in order for c in conv_n[u].values()] or [1])
    norm = Normalize(vmin=1, vmax=max(args.vmax or maxc, 2))

    fig_h = min(max(3.0, 0.16 * n + 2.0), 30)
    fig, ax = plt.subplots(figsize=(args.width, fig_h))

    # SNP mask positions as faint vertical lines
    for p in snps:
        if x0 <= p + 1 <= x1:
            ax.axvline(p + 1, color="orange", lw=0.6, alpha=0.5, zorder=0)

    for i, u in enumerate(order):
        y = i + 1
        ax.axhline(y, color="0.85", lw=0.4, ls=":", zorder=0)
        absent = [p + 1 for p in covered[u] if conv_n[u][p] == 0]
        ax.scatter(absent, [y] * len(absent), s=args.dot_size, facecolors="0.75",
                   edgecolors="0.55", linewidths=0.3, zorder=1)
        pres = [(p, c) for p, c in conv_n[u].items() if p not in snps]
        if pres:
            ax.scatter([p + 1 for p, _ in pres], [y] * len(pres),
                       c=[c for _, c in pres], cmap=cmap, norm=norm, marker="x",
                       s=args.dot_size * 2.2, linewidths=1.3, zorder=3)
        msk = [p + 1 for p in conv_n[u] if p in snps]
        if msk:
            ax.scatter(msk, [y] * len(msk), marker="s", s=args.dot_size * 2,
                       facecolors="none", edgecolors="darkorange",
                       linewidths=1.0, zorder=2)

    ax.set_xlim(x0, x1)
    ax.set_ylim(n + 0.8, 0.2)
    ax.set_yticks([1, n] if n > 1 else [1])
    ax.set_ylabel(f"{n} transcripts ({len([u for u in order if u in lab])} labeled, "
                  f"{len([u for u in order if u in unl])} unlabeled)")
    ax.set_xlabel(f"{chrom} position")
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.grid(axis="x", color="0.9", lw=0.8)
    ax.set_title(f"{gene} ({gene_strand} strand) — cell {cell}")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.03,
                      pad=0.08 if n > 20 else 0.18, aspect=40)
    cb.set_label("Reads supporting T-to-C per uridine site")

    handles = [
        Line2D([], [], marker="o", ls="", mfc="0.75", mec="0.55", label="T-to-C absent"),
        Line2D([], [], marker="x", ls="", color="#0000ff", mew=1.5, label="T-to-C present"),
    ]
    if snps:
        handles += [
            Line2D([], [], marker="s", ls="", mfc="none", mec="darkorange",
                   label="T-to-C at masked SNP (ignored)"),
            Line2D([], [], color="orange", alpha=0.6, label="SNP mask position"),
        ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"\nSaved plot: {args.out}  ({len(lab)} labeled / {len(unl)} unlabeled UMIs)")
    return order


def write_table(path, order, covered, conv_n, snps, chrom):
    with open(path, "w") as f:
        f.write("umi\tchrom\tpos_1based\tTC_reads\tin_snp_mask\n")
        for u in order:
            for p in sorted(covered[u]):
                f.write(f"{u}\t{chrom}\t{p + 1}\t{conv_n[u][p]}\t{int(p in snps)}\n")
    print(f"Saved table: {path}")


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="pre-consensus, coordinate-sorted, indexed BAM")
    ap.add_argument("--fasta", required=True, help="genome FASTA (indexed, .fai)")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--gene", required=True, help="gene_name, e.g. FOS")
    ap.add_argument("--cell", help="cell barcode (without -1). Default: cell with most UMIs")
    ap.add_argument("--snp-csv", help="SNP mask CSV (contig,genome_i[,conversion])")
    ap.add_argument("--snp-coord", choices=["1based", "0based"], default="1based",
                    help="coordinate system of genome_i in --snp-csv (VCF POS = 1based)")
    ap.add_argument("--strand", choices=["forward", "reverse", "unstranded"],
                    default="forward", help="library strandedness (10x 3' = forward)")
    ap.add_argument("--cb-tag", default="CB")
    ap.add_argument("--ub-tag", default="UB")
    ap.add_argument("--min-bq", type=int, default=20, help="min base quality")
    ap.add_argument("--min-mapq", type=int, default=255,
                    help="min MAPQ (STAR unique = 255)")
    ap.add_argument("--min-conv-reads", type=int, default=1,
                    help="reads supporting a T-to-C for a UMI to count as labeled")
    ap.add_argument("--max-umis", type=int, default=200)
    ap.add_argument("--xlim", type=int, nargs=2, help="1-based x range")
    ap.add_argument("--vmax", type=float, help="colour scale max (default: data max)")
    ap.add_argument("--qc-frac", type=float, default=0.5,
                    help="flag sites converted in >= this fraction of covering UMIs")
    ap.add_argument("--qc-min-umis", type=int, default=5)
    ap.add_argument("--width", type=float, default=7.5)
    ap.add_argument("--dot-size", type=float, default=9)
    ap.add_argument("--out", default="tc_reads.pdf")
    ap.add_argument("--table", help="optional TSV of per-UMI per-uridine calls")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.fasta)
    chrom, start, end, gstrand, gid = find_gene(args.gtf, args.gene, set(bam.references))
    print(f"{args.gene} ({gid}): {chrom}:{start + 1}-{end} ({gstrand})")

    cell = norm_bc(args.cell) if args.cell else rank_cells(
        bam, chrom, start, end, gstrand, args)
    print(f"Cell: {cell}")

    snps = load_snps(args.snp_csv, chrom, args.snp_coord) if args.snp_csv else set()
    snps = {p for p in snps if start <= p < end}
    if args.snp_csv:
        print(f"SNP mask sites inside {args.gene}: {len(snps)}")

    covered, conv_n, n_reads = collect(bam, fasta, chrom, start, end, gstrand, cell, args)
    if n_reads and np.median(list(n_reads.values())) <= 1:
        print("NOTE: median 1 read per UMI -- if this is a dynast consensus BAM, "
              "use the pre-consensus BAM to get meaningful coverage colours.")

    snp_qc(covered, conv_n, snps, chrom, args.qc_frac, args.qc_min_umis)
    order = plot(covered, conv_n, snps, chrom, args.gene, gstrand, cell, args)
    if args.table:
        write_table(args.table, order, covered, conv_n, snps, chrom)


if __name__ == "__main__":
    main()
