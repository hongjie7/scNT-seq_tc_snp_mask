#!/usr/bin/env python3
"""
plot_tc_panel.py -- scNT-seq style T-to-C plots of single transcripts (UMIs) in one cell.

Every row is one transcript (UMI). Grey circle = covered uridine without T-to-C;
X = uridine with T-to-C, coloured by how many reads of that UMI carry it.
Masked SNP positions (--snp-csv, 0-based) are not counted. The cell barcode is
added to the output file name and printed on the figure.

INPUT  --bam = pre-consensus BAM (split_bams/<donor>.bam, input of dynast consensus)
       --compare-consensus = consensus.bam (output of dynast consensus)
       --fasta, --gtf, --genes are always required.

PLOT TYPES (choose with the options in the "Plot type" group)
 1 basic         one panel per gene, same cell (fast vs stable gene, Fig. 1c)
                   --genes MYC EEF1A1 --labels "fast" "stable"   [--top-cells 3]
                   or genes by NTR: --adata <masked adata.h5ad> [--high-ntr-genes 1 --low-ntr-genes 1]
                   every read inside its row: --reads-in-rows (basic_allreads)
 2 before-after  one gene: no mask | (negative controls) | corrected mask;
                 zooms on the masked SNP that removes most T>C
                   --genes GSTP1 --before-after [--only-masked]
                   [--control-shift 1] [--control-mask other_donor.csv]
                   [--top-cells 3]   (one figure per cell at the same SNP)
                 several genes -> one figure per gene; or pick genes automatically:
                   --mask-genes 3 --adata <masked adata.h5ad> --adata-unmasked <TC adata.h5ad>
                   [--min-removed 5] [--min-removed-pct 20] [--mask-rank-by removed|pct|umis|ntr]
                   [--list-mask-genes 20]   (table only)
 3 scnt_layout   scNT-seq (Qiu 2020) ACTG1 figure layout: gene bar, transcript
                 rows, highlighted transcript and an inset with all reads of its UMI
                   --genes ACTG1 --scnt-layout [--utr3] [--umi-reads UMI] [--label-umis]
 4 scnt_layout, pre vs consensus
                   --scnt-layout --compare-consensus consensus.bam --compare-layout inset
                     (inset: reads + consensus record of the UMI)
                   --scnt-layout --compare-consensus consensus.bam --compare-layout side
                     (whole figure twice: pre-consensus | consensus.bam)
 5 panels, pre vs consensus
                   --compare-consensus consensus.bam --compare-panels 2x2|masked|nomask
                   [--highlight-umis U1 U2] [--reads-in-rows] [--label-umis]
 6 one UMI, pre vs consensus (reads on top, consensus record below)
                   --compare-consensus consensus.bam --umi-compare [--umi-reads UMI]
 7 reads of one UMI (pre-consensus BAM)
                   --umi-reads UMI [--umi-only]   (panel below the standard plot)
 8 list cells    --list-cells 5   (top cells per gene, no plot)

Which BAM: report counts from consensus.bam (--min-bq 28 = dynast's cutoff);
use the pre-consensus BAM for read-level views and coverage colours.

EXAMPLES
  bam=split_bams/DONOR.bam; mask=01_mask/DONOR_snps_0based.csv
  common="--fasta 00_input/genome.fa --gtf 00_input/annotation.gtf --snp-csv $mask"

  python plot_tc_panel.py --bam $bam $common --genes MYC EEF1A1 --vmax 4 --out fast_vs_stable.pdf
  python plot_tc_panel.py --bam 00_input/consensus.bam $common --genes GSTP1 \\
      --before-after --min-bq 28 --control-shift 1 --window 50 --out GSTP1_mask.pdf
  python plot_tc_panel.py --bam $bam $common --genes ACTG1 --utr3 --scnt-layout \\
      --label-umis --sample-label "Day 0" --out ACTG1_paper.pdf
  python plot_tc_panel.py --bam $bam $common --genes SOX2 --scnt-layout --label-umis \\
      --compare-consensus 00_input/consensus.bam --compare-layout side --out SOX2_side.pdf
  python plot_tc_panel.py --bam $bam $common --genes ACTB \\
      --compare-consensus 00_input/consensus.bam --compare-panels masked \\
      --label-umis --max-umis 100000 --reads-in-rows --row-height 0.3 --out ACTB_cmp.pdf
  python plot_tc_panel.py --bam $bam $common --genes ACTB --umi-compare \\
      --compare-consensus 00_input/consensus.bam --out ACTB_umi.pdf
  python plot_tc_panel.py --bam $bam $common --genes SOX2 ACTB --list-cells 5

Requires plot_tc_reads.py in the same folder (reuses its read parsing).
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_tc_reads import collect, load_snps, norm_bc, iter_gene_reads  # noqa: E402

def out_with_cell(path, cells):
    """Insert the cell barcode(s) into the output file name: plot.pdf -> plot_<cell>.pdf.
    Not repeated if the name already contains the barcode."""
    cells = [c for c in dict.fromkeys(cells) if c]
    root, ext = os.path.splitext(path)
    tag = "_".join(c for c in cells if c not in os.path.basename(root))
    return f"{root}_{tag}{ext}" if tag else path


PLOT_NAME = "basic"      # set in main() from the chosen plot type


def plot_name(a):
    """Short name of the plot type, used in file names and on the figure."""
    if a.umi_compare:
        return "umi_compare"
    if a.paper and a.compare_layout == "mask":
        return "inset_mask"
    if a.paper and a.compare_consensus:
        return ("scnt_layout_consensus_side" if a.compare_layout == "side"
                else "scnt_layout_consensus_inset")
    if a.paper:
        return "scnt_layout"
    if a.compare_consensus:
        return "consensus_" + a.compare_panels
    if a.before_after:
        return "mask_only" if a.only_masked else "mask_before_after"
    if a.umi_reads and a.umi_only:
        return "umi_reads"
    return "basic_allreads" if a.reads_in_rows else "basic"


def name_tags(a):
    """Settings written into the output file name: group-yes/no and mask-yes/no.
    (before/after and consensus panel plots already say the mask in their plot name)"""
    tags = ["group-yes" if a.order == "labeled_first" else "group-no"]
    if not (PLOT_NAME.startswith(("mask_", "consensus_")) or PLOT_NAME == "inset_mask"):
        tags.append("mask-yes" if a.snp_csv else "mask-no")
    return tags


def add_tags(path, tags):
    """For a user --out name: insert the tags before the extension unless the name already
    contains a group-/mask- tag."""
    root, ext = os.path.splitext(path)
    base = os.path.basename(root)
    new = [t for t in tags if t.split("-")[0] + "-" not in base]
    return f"{root}{''.join('_' + t for t in new)}{ext}"


def label_cells(fig, cells):
    """Write the plot type and cell barcode(s) in the top-left corner of the figure."""
    cells = [c for c in dict.fromkeys(cells) if c]
    if cells:
        grp = "yes" if ORDER == "labeled_first" else "no"
        fig.text(0.01, 0.995, f"{PLOT_NAME} plot  |  group by label: {grp}  |  Cell: "
                 + ", ".join(cells), ha="left",
                 va="top", fontsize=8, color="0.3")


CMAP = LinearSegmentedColormap.from_list("bpr", ["#1f1fff", "#8a2be2", "#e8302a"])


def find_genes(gtf, names, bam_contigs):
    """One pass over the GTF -> {name: (chrom, start0, end, strand)}."""
    want = set(names)
    hits = defaultdict(list)
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as f:
        for line in f:
            if line.startswith("#") or 'gene_name "' not in line:
                continue
            p = line.split("\t", 8)
            if p[2] != "gene":
                continue
            name = p[8].split('gene_name "')[1].split('"')[0]
            if name in want:
                hits[name].append((p[0], int(p[3]) - 1, int(p[4]), p[6]))
    out = {}
    for n in names:
        h = sorted(hits.get(n, []), key=lambda x: (x[0] not in bam_contigs, "_" in x[0]))
        if not h or h[0][0] not in bam_contigs:
            sys.exit(f"ERROR: gene '{n}' not found in GTF on a BAM contig "
                     "(gene names are case-sensitive: human FOS, mouse Fos)")
        out[n] = h[0]
    return out


def utr3_span(gtf, name):
    """1-based (start, end) of the 3'UTR of the gene's canonical transcript
    (tag Ensembl_canonical or MANE_Select; else the transcript with the longest 3'UTR).
    Works with GENCODE ("UTR" + stop_codon), Ensembl ("three_prime_utr") and
    "three_prime_UTR" feature names; falls back to exon minus CDS if no UTR lines."""
    opener = gzip.open if gtf.endswith(".gz") else open
    key = f'gene_name "{name}";'
    utr3, utr, stop, canon = defaultdict(list), defaultdict(list), {}, set()
    exons, cds, strand_of, ftypes = defaultdict(list), defaultdict(list), {}, set()
    with opener(gtf, "rt") as f:
        for line in f:
            if key not in line or line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if "_" in p[0] or 'transcript_id "' not in p[8]:
                continue
            tid = p[8].split('transcript_id "')[1].split('"')[0]
            ft, st, a0, b0 = p[2].lower(), p[6], int(p[3]), int(p[4])
            ftypes.add(p[2]); strand_of[tid] = st
            if ft == "transcript":
                if '"Ensembl_canonical"' in p[8] or '"MANE_Select"' in p[8]:
                    canon.add(tid)
            elif ft in ("three_prime_utr", "3utr", "3'utr"):
                utr3[tid].append((a0, b0))
            elif ft == "utr":
                utr[tid].append((a0, b0))
            elif ft == "stop_codon":
                o = stop.get(tid)
                stop[tid] = (min(a0, o[0]) if o else a0, max(b0, o[1]) if o else b0)
            elif ft == "exon":
                exons[tid].append((a0, b0))
            elif ft == "cds":
                cds[tid].append((a0, b0))
    spans = {}
    for tid in set(utr3) | set(utr) | set(exons):
        st = strand_of[tid]
        three = list(utr3.get(tid, []))
        if not three and tid in utr and (tid in stop or tid in cds):
            # GENCODE: generic "UTR" -> keep pieces 3' of the stop codon / CDS
            if tid in stop:
                s0, e0 = stop[tid]
            else:
                s0, e0 = min(x for x, _ in cds[tid]), max(y for _, y in cds[tid])
            three = [(x, y) for x, y in utr[tid] if (x > e0 if st == "+" else y < s0)]
        if not three and tid in exons and tid in cds:
            # fallback: exon parts beyond the CDS (stop codon included in CDS or not)
            c0, c1 = min(x for x, _ in cds[tid]), max(y for _, y in cds[tid])
            if tid in stop:
                c0, c1 = min(c0, stop[tid][0]), max(c1, stop[tid][1])
            for x, y in exons[tid]:
                if st == "+" and y > c1:
                    three.append((max(x, c1 + 1), y))
                elif st == "-" and x < c0:
                    three.append((x, min(y, c0 - 1)))
        if three:
            spans[tid] = (min(x for x, _ in three), max(y for _, y in three))
    if not spans and exons:
        # no UTR/CDS information at all: use the 3'-most exon of each transcript
        for tid, ex in exons.items():
            x, y = max(ex) if strand_of[tid] == "+" else min(ex)
            spans[tid] = (x, y)
        print(f"NOTE: the GTF has no UTR/CDS lines for {name}; showing the 3'-most exon instead")
    if not spans:
        sys.exit(f"ERROR: no 3'UTR found for {name} in the GTF (feature types seen: "
                 f"{', '.join(sorted(ftypes)) or 'none'}). Use --xlim instead.")
    pick = [t for t in spans if t in canon] or sorted(spans, key=lambda t: spans[t][0] - spans[t][1])
    tid = pick[0]
    print(f"3'UTR from transcript {tid}" + (" (canonical)" if tid in canon else
                                           " (no canonical tag; longest 3'UTR)"))
    return spans[tid]


def umis_per_cell(bam, loc, args):
    chrom, start, end, strand = loc
    d = defaultdict(set)
    for r in iter_gene_reads(bam, chrom, start, end, strand, args.strand,
                             args.cb_tag, args.ub_tag, args.min_mapq):
        d[norm_bc(r.get_tag(args.cb_tag))].add(r.get_tag(args.ub_tag))
    return {c: len(u) for c, u in d.items()}


def report_genes(per_gene, genes, top=5):
    for g in genes:
        d = per_gene[g]
        tot = sum(d.values())
        best = sorted(d.items(), key=lambda kv: -kv[1])[:top]
        print(f"  {g}: {len(d):,} cells with >=1 UMI, {tot:,} UMIs in total; top cells: "
              + (", ".join(f"{c}={n}" for c, n in best) or "none"))


def pick_cell(per_gene, genes, top=10):
    cells = set.intersection(*[set(per_gene[g]) for g in genes])
    if not cells:
        print("No cell has UMIs in all listed genes:")
        report_genes(per_gene, genes)
        sys.exit("ERROR: use --cell-mode each (best cell per gene), other genes, or another sample")
    ranked = sorted(cells, key=lambda c: (min(per_gene[g][c] for g in genes),
                                          sum(per_gene[g][c] for g in genes)), reverse=True)
    print("Top cells (UMIs per gene):")
    for c in ranked[:top]:
        print("  " + c + "\t" + "\t".join(f"{g}={per_gene[g][c]}" for g in genes))
    return ranked[0]


def scan_removed(bam, fa, loc, snps, args):
    """For a --before-after plot: find the masked SNP sites in the gene where the
    mask removes the most T>C, and per site the cells with the most affected UMIs.
    Returns {site0: {"reads": n, "umis": n, "cells": {cell: [umis_conv, umis_cov]}}}."""
    chrom, start, end, strand = loc
    conv_b = "C" if strand == "+" else "G"
    conv_u = defaultdict(set)    # (site, cell) -> UMIs with T>C at site
    cov_u = defaultdict(set)     # (site, cell) -> UMIs covering site
    reads = defaultdict(int)
    for r in iter_gene_reads(bam, chrom, start, end, strand, args.strand,
                             args.cb_tag, args.ub_tag, args.min_mapq):
        if r.reference_end is None:
            continue
        hit = [p for p in snps if r.reference_start <= p < r.reference_end] if len(snps) < 50 \
            else None
        if hit is not None and not hit:
            continue
        cell, umi = norm_bc(r.get_tag(args.cb_tag)), r.get_tag(args.ub_tag)
        seq, qual = r.query_sequence, r.query_qualities
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp not in snps:
                continue
            if qual is not None and qual[qp] < args.min_bq:
                continue
            cov_u[(rp, cell)].add(umi)
            if seq[qp] == conv_b:
                conv_u[(rp, cell)].add(umi)
                reads[rp] += 1
    out = {}
    for (site, cell), us in conv_u.items():
        d = out.setdefault(site, {"reads": reads[site], "umis": 0, "cells": {}})
        d["umis"] += len(us)
        d["cells"][cell] = [len(us), len(cov_u[(site, cell)])]
    return out


SUBSAMPLE = "labeled_first"   # set from --subsample in main()
ORDER = "labeled_first"       # set from --order in main()


def collect_umi_reads(bam, fa, loc, cell, umi, args, snps):
    """Reads of one UMI in one cell, each as (aligned blocks, uridines, T>C sites)."""
    chrom, start, end, strand = loc
    ref_u, conv_b = ("T", "C") if strand == "+" else ("A", "G")
    out = []
    for r in iter_gene_reads(bam, chrom, start, end, strand, args.strand,
                             args.cb_tag, args.ub_tag, args.min_mapq):
        if norm_bc(r.get_tag(args.cb_tag)) != cell or r.get_tag(args.ub_tag) != umi:
            continue
        ref = fa.fetch(chrom, r.reference_start, r.reference_end).upper()
        seq, qual = r.query_sequence, r.query_qualities
        uri, tc = set(), set()
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if ref[rp - r.reference_start] != ref_u or seq[qp] == "N":
                continue
            if qual is not None and qual[qp] < args.min_bq:
                continue
            uri.add(rp)
            if seq[qp] == conv_b and rp not in snps:
                tc.add(rp)
        out.append((r.get_blocks(), uri, tc))
    out.sort(key=lambda x: (x[0][0][0], x[0][-1][1]) if x[0] else (0, 0))
    return out


def collect_cell_reads(bam, fa, loc, cell, args):
    """All reads of one cell for a gene, grouped by UMI: {umi: [(start, end, uridines, tc)]}.
    For a dynast consensus BAM, original single reads duplicated by a RN:i:1 copy are dropped."""
    chrom, start, end, strand = loc
    ref_u, conv_b = ("T", "C") if strand == "+" else ("A", "G")
    by = defaultdict(list)
    for r in iter_gene_reads(bam, chrom, start, end, strand, args.strand,
                             args.cb_tag, args.ub_tag, args.min_mapq):
        if norm_bc(r.get_tag(args.cb_tag)) == cell:
            by[r.get_tag(args.ub_tag)].append(r)
    out = {}
    for umi, rs in by.items():
        cons = [r for r in rs if r.has_tag("RN")]
        rows = []
        for r in (cons or rs):
            ref = fa.fetch(chrom, r.reference_start, r.reference_end).upper()
            seq, qual = r.query_sequence, r.query_qualities
            uri, tc = set(), set()
            for qp, rp in r.get_aligned_pairs(matches_only=True):
                if rp < start or rp >= end or ref[rp - r.reference_start] != ref_u:
                    continue
                if seq[qp] == "N" or (qual is not None and qual[qp] < args.min_bq):
                    continue
                uri.add(rp)
                if seq[qp] == conv_b:
                    tc.add(rp)
            if uri:
                rows.append((r.reference_start, r.reference_end, uri, tc))
        rows.sort()
        out[umi] = rows
    return out


def draw_umi_reads(ax, reads, chrom, umi, args, xlim=None):
    """Paper-style inset: one line per read of a single UMI; dotted line = gap (splice / N)."""
    for i, (blocks, uri, tc) in enumerate(reads, 1):
        if blocks:
            ax.plot([blocks[0][0] + 1, blocks[-1][1]], [i, i], ls=(0, (1, 1.5)),
                    color="0.6", lw=0.8, zorder=0)
        ab = [p + 1 for p in uri if p not in tc]
        ax.scatter(ab, [i] * len(ab), s=args.dot_size, facecolors="0.86", edgecolors="0.55",
                   linewidths=0.35, zorder=1)
        if tc:
            ax.scatter([p + 1 for p in tc], [i] * len(tc), marker="x", color="#1f1fff",
                       s=args.dot_size * 2.6, linewidths=1.4, zorder=3)
    n = len(reads)
    ax.set_ylim(n + 0.7, 0.3)
    ax.set_yticks(range(1, n + 1))
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=3)
    if xlim:
        ax.set_xlim(*xlim)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5, integer=True))
    ax.grid(axis="x", color="0.85", lw=1.0, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel(f"{chrom} position", fontsize=9, labelpad=2)
    ax.set_title(f"{n} reads sharing the same UMI index ({umi})", fontsize=10)


PAPER_COLORS = ["#1f1fff", "#9b2bb5", "#cc2b5e", "#e8302a"]     # coverage 1, 2, 3, >=4


def plot_paper(a, g, loc, cell, cov, conv, snps, cell_reads, region_label, cons_reads=None,
               fig=None, rep_override=None, order_override=None, xlim_override=None,
               panel_title=None):
    """Layout of scNT-seq Fig. (ACTG1 3'UTR): gene bar with strand arrows, one row per
    transcript (UMI), discrete coverage colours 1-4, a highlighted representative transcript
    and an inset with all reads sharing its UMI."""
    from matplotlib.patches import Rectangle
    chrom, gstart, gend, strand = loc

    order, lab, unl = order_umis(cov, conv, snps, a.min_conv_reads, a.max_umis)
    if order_override is not None:
        order = list(order_override)
    n = len(order)
    if n == 0:
        sys.exit(f"ERROR: no transcripts of {g} in cell {cell} in this window")
    allpos = np.array([p for u in order for p in cov[u]] or [loc[1]]) + 1
    if xlim_override:
        x0, x1 = xlim_override
    elif a.xlim:
        x0, x1 = a.xlim
    else:
        lo, hi = np.percentile(allpos, [a.range_pct, 100 - a.range_pct])
        pad = max(10, int(0.03 * (hi - lo)))
        x0, x1 = int(lo) - pad, int(hi) + pad

    # representative transcript: labeled, most reads
    if rep_override:
        rep = rep_override
    elif a.umi_reads and a.umi_reads != "auto":
        rep = a.umi_reads
    else:
        cand = [u for u in order if cell_reads.get(u)]
        rep = max(cand, key=lambda u: (u in lab, len(cell_reads[u]))) if cand else None
    rep_reads = cell_reads.get(rep, []) if rep else []

    col = lambda c: PAPER_COLORS[min(c, 4) - 1]
    row_h = a.row_height or 0.13
    main_h = max(2.5, row_h * n + 0.6)
    cons_rep = (cons_reads or {}).get(rep, []) if rep else []
    n_inset = len(rep_reads) + (len(cons_rep) + 1 if cons_reads is not None else 0)
    inset_h = max(2.6, 0.22 * n_inset + 0.9)
    standalone = fig is None
    if standalone:
        fig = plt.figure(figsize=(a.width, 0.9 + main_h + inset_h + 0.6))
    else:
        main_h, inset_h = a._paper_h       # same sizes in both side-by-side figures
    root = fig.figure if hasattr(fig, "figure") and fig.figure is not None else fig
    # hspace is relative to the mean panel height: keep the gaps at ~1 inch
    gs = fig.add_gridspec(3, 2, height_ratios=[0.55, main_h, inset_h], width_ratios=[1, 1.15],
                          hspace=(1.25 if (cons_reads is not None or not standalone) else 1.0)
                          / np.mean([0.55, main_h, inset_h]),
                          wspace=0.15)

    # gene bar with strand arrows
    bax = fig.add_subplot(gs[0, :]); bax.set_xlim(0, 1); bax.set_ylim(0, 1); bax.axis("off")
    bax.text(0.5, 0.92, f"$\\it{{{g}}}$ {region_label} ({chrom}: {x0:,}-{x1:,})",
             ha="center", va="top", fontsize=11)
    bax.add_patch(Rectangle((0.08, 0.05), 0.84, 0.35, color="0.3"))
    bax.text(0.5, 0.22, ("<" if strand == "-" else ">") * 4, ha="center", va="center",
             color="#9db4ff", fontsize=13, fontweight="bold")

    # main panel
    ax = fig.add_subplot(gs[1, :])
    rep_row = None
    for i, u in enumerate(order, 1):
        ax.axhline(i, color="0.9", lw=0.4, ls=(0, (1, 2)), zorder=0)
        sites = {p: c for p, c in conv[u].items() if p not in snps and c >= a.min_conv_reads}
        msk = [p + 1 for p, c in conv[u].items() if p in snps] if a.show_masked else []
        ab = [p + 1 for p in cov[u] if p not in sites and p + 1 not in msk]
        if msk:
            ax.scatter(msk, [i] * len(msk), marker="s", s=a.dot_size * 2.2, facecolors="none",
                       edgecolors="darkorange", linewidths=1.0, zorder=2)
        ax.scatter(ab, [i] * len(ab), s=a.dot_size, facecolors="0.88", edgecolors="0.55",
                   linewidths=0.4, zorder=1)
        if sites:
            ax.scatter([p + 1 for p in sites], [i] * len(sites),
                       c=[col(c) for c in sites.values()], marker="x",
                       s=a.dot_size * 2.6, linewidths=1.5, zorder=3)
        if u == rep:
            rep_row = i
    ax.set_xlim(x0, x1); ax.set_ylim(n + 0.7, 0.3)
    if a.label_umis:
        # UMI sequence of every row; the representative transcript in bold red
        ax.set_yticks(range(1, n + 1))
        ax.set_yticklabels(order, fontsize=a.umi_fontsize, family="monospace")
        for t, u in zip(ax.get_yticklabels(), order):
            if u == rep:
                t.set_color("#e8302a"); t.set_fontweight("bold")
    else:
        ticks = [1] + [t for t in range(5, n, 5)] + [n] if n > 6 else list(range(1, n + 1))
        ax.set_yticks(ticks); ax.set_yticklabels([str(t) if t in (1, n) else "" for t in ticks])
    ax.set_xlabel(f"{chrom} position", fontsize=9, labelpad=2)
    ax.set_ylabel(f"{n} transcripts ({len(lab & set(order))} labeled, "
                  f"{len(unl & set(order))} unlabeled)", fontsize=10)
    if n < len(lab) + len(unl):
        ax.set_ylabel(f"{n} of {len(lab) + len(unl)} transcripts\n"
                      f"(all: {len(lab)} labeled, {len(unl)} unlabeled)", fontsize=10)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5, integer=True))
    ax.grid(axis="x", color="0.85", lw=1.0, zorder=0); ax.set_axisbelow(True)
    if panel_title or a.sample_label:
        ax.set_title(panel_title or a.sample_label, fontsize=11)
    for sp in ax.spines.values():
        sp.set_linewidth(1.2)

    # legend (bottom left)
    lax = fig.add_subplot(gs[2, 0]); lax.axis("off")
    lax.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc="0.88", mec="0.55", ms=9, label="Absent"),
        Line2D([], [], marker="x", ls="", color=PAPER_COLORS[0], mew=2.2, ms=9, label="Present")] + ([Line2D([], [], marker="s", ls="", mfc="none", mec="darkorange", ms=9, label="Masked SNP (not counted)")] if a.show_masked else []),
        title="T-to-C substitution at\na uridine site:", loc="upper left",
        bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=10, title_fontsize=10)
    lax.legend_.set_title("T-to-C substitution at\na uridine site:")
    leg1 = lax.legend_
    lax.add_artist(leg1)
    lax.legend(handles=[Line2D([], [], marker="s", ls="", color=c, ms=10,
                               label=("≥4" if k == 4 else str(k)))
                        for k, c in enumerate(PAPER_COLORS, 1)],
               title="Coverage of T-to-C\nsubstitutions per uridine site:", loc="lower left",
               bbox_to_anchor=(0.0, 0.0), frameon=False, fontsize=10, title_fontsize=10,
               ncol=4, handletextpad=0.2, columnspacing=0.8)

    # inset: reads of the representative transcript
    iax = fig.add_subplot(gs[2, 1])
    if rep_reads:
        k = len(rep_reads)
        rows = [("read", r) for r in rep_reads]
        if cons_reads is not None:
            rows.append(("sep", None))
            rows += [("cons", r) for r in cons_rep]
        ylabels = []
        for j, (kind, rr) in enumerate(rows, 1):
            if kind == "sep":
                iax.axhline(j, color="0.4", lw=0.8)
                ylabels.append("")
                continue
            r0, r1, uri, tc = rr
            colr = "#d97706" if kind == "cons" else "0.45"
            iax.plot([r0 + 1, r1], [j, j], ls=(0, (1.5, 1.5)), color=colr, lw=0.9, zorder=0)
            ab = [p + 1 for p in uri if p not in tc or (p in snps and not a.show_masked)]
            iax.scatter(ab, [j] * len(ab), s=a.dot_size, facecolors="0.88", edgecolors="0.45",
                        linewidths=0.4, zorder=1)
            if a.show_masked:
                ms = [p + 1 for p in tc if p in snps]
                if ms:
                    iax.scatter(ms, [j] * len(ms), marker="s", s=a.dot_size * 2.2,
                                facecolors="none", edgecolors="darkorange", linewidths=1.0,
                                zorder=2)
            xs = [p + 1 for p in tc if p not in snps]
            if xs:
                iax.scatter(xs, [j] * len(xs), marker="x", color=PAPER_COLORS[0],
                            s=a.dot_size * 2.6, linewidths=1.5, zorder=3)
            ylabels.append(f"read {j}" if kind == "read" else "consensus")
        if cons_reads is not None and not cons_rep:
            iax.text(0.02, (k + 1.5) / (len(rows) + 1), "UMI not in consensus.bam (window/filters)",
                     transform=iax.transAxes, fontsize=7, color="0.4")
        allr = rep_reads + list(cons_rep)
        lo_ = min(r[0] for r in allr) + 1; hi_ = max(r[1] for r in allr)
        pad = max(8, int(0.08 * (hi_ - lo_)))
        iax.set_xlim(lo_ - pad, hi_ + pad); iax.set_ylim(len(rows) + 0.7, 0.3)
        iax.set_yticks(range(1, len(rows) + 1))
        iax.set_yticklabels(ylabels if cons_reads is not None else [], fontsize=7)
        for t, lab_ in zip(iax.get_yticklabels(), ylabels):
            if lab_ == "consensus":
                t.set_color("#d97706"); t.set_fontweight("bold")
        iax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
        iax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3, integer=True))
        iax.grid(axis="x", color="0.85", lw=1.0, zorder=0); iax.set_axisbelow(True)
        iax.set_title(f"{k} read{'s' if k != 1 else ''} sharing the same UMI index\n"
                      f"for a representative transcript"
                      , fontsize=10)
        iax.text(0.5, -0.12 if len(rows) > 6 else -0.25, f"UMI {rep}", transform=iax.transAxes,
                 ha="center", va="top", fontsize=9, color="#e8302a", family="monospace",
                 fontweight="bold")
        for sp in iax.spines.values():
            sp.set_linewidth(1.2)
        if k <= 1:
            print("NOTE: the representative transcript has 1 read -- use the pre-consensus BAM")
        # highlight the row and connect it to the inset
        if rep_row:
            ax.add_patch(Rectangle((x0, rep_row - 0.5), x1 - x0, 1.0, fill=False, ls="--",
                                   ec="#e8302a", lw=1.2, zorder=4, clip_on=False))
            # bracket arrow: row (right end) -> right margin -> down -> into the inset
            to_fig = root.transFigure.inverted()
            ax_box = ax.get_window_extent().transformed(to_fig)
            in_box = iax.get_window_extent().transformed(to_fig)
            y_row = to_fig.transform(ax.transData.transform((x1, rep_row)))[1]
            y_in = in_box.y0 + 0.5 * in_box.height
            x_out = max(ax_box.x1, in_box.x1) + (0.035 if standalone else 0.015)
            root.add_artist(Line2D([ax_box.x1, x_out, x_out], [y_row, y_row, y_in],
                                   transform=root.transFigure, color="black", lw=1.2))
            root.add_artist(matplotlib.patches.FancyArrowPatch(
                (x_out, y_in), (in_box.x1 + 0.003, y_in), transform=root.transFigure,
                arrowstyle="-|>", mutation_scale=12, color="black", lw=1.2))
        print(f"Representative transcript: UMI {rep} ({k} reads), row {rep_row}")
    else:
        iax.axis("off")
    print(f"{g}{' (' + panel_title + ')' if panel_title else ''}: {n} transcripts shown, "
          f"{len(lab)} labeled / {len(unl)} unlabeled in cell {cell}")
    if not standalone:
        return rep
    label_cells(fig, [cell])
    a.out = out_with_cell(a.out, [cell])
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    print(f"Saved: {a.out}")


def order_umis(covered, conv_n, snps, min_conv_reads, max_umis):
    """Rows: labeled transcripts first, then unlabeled; within each group sorted by the
    first covered position (left to right), like the scNT-seq figure.
    If there are more than max_umis transcripts, which ones are kept depends on --subsample:
      labeled_first : all labeled (up to max_umis), unlabeled filled up by RANDOM choice
      proportional  : random sample keeping the labeled fraction of the cell
      first         : the first max_umis after sorting (old behaviour; biased to the left)"""
    import random

    def labeled(u):
        return any(n >= min_conv_reads and p not in snps for p, n in conv_n[u].items())
    umis = sorted(u for u in covered if covered[u])
    lab = [u for u in umis if labeled(u)]
    unl = [u for u in umis if not labeled(u)]
    start = lambda u: min(covered[u])
    if len(umis) > max_umis:
        rng = random.Random(0)                     # reproducible choice
        if SUBSAMPLE == "proportional":
            k = round(max_umis * len(lab) / len(umis))
            lab_k, unl_k = rng.sample(lab, k), rng.sample(unl, min(len(unl), max_umis - k))
        elif SUBSAMPLE == "labeled_first":
            lab_k = lab if len(lab) <= max_umis else rng.sample(lab, max_umis)
            unl_k = rng.sample(unl, min(len(unl), max_umis - len(lab_k)))
        else:
            lab_k, unl_k = sorted(lab, key=start)[:max_umis], None
            unl_k = sorted(unl, key=start)[:max(0, max_umis - len(lab_k))]
    else:
        lab_k, unl_k = lab, unl
    if ORDER == "position":       # all rows by first covered position, labeled mixed in
        order = sorted(lab_k + unl_k, key=lambda u: (start(u), max(covered[u])))
    else:                         # labeled block on top, each block by position
        order = sorted(lab_k, key=start) + sorted(unl_k, key=start)
    return order, set(lab), set(unl)


def draw_panel(ax, gene, label, loc, covered, conv_n, snps, args, norm,
               order=None, title=None, mark=None):
    chrom = loc[0]
    order0, lab, unl = order_umis(covered, conv_n, snps, args.min_conv_reads, args.max_umis)
    if order is None:
        order = order0
    n = len(order)
    if n == 0:
        ax.text(0.5, 0.5, f"{gene}: no reads in this cell", ha="center", transform=ax.transAxes)
        return 0, 0
    allpos = np.array([p for u in order for p in covered[u]]) + 1
    if args.full_range:
        lo, hi = allpos.min(), allpos.max()
    else:   # ignore a few far-off reads (e.g. upstream exons) that squeeze the main signal
        lo, hi = np.percentile(allpos, [args.range_pct, 100 - args.range_pct])
    pad = max(15, int(0.03 * (hi - lo)))
    x0, x1 = int(lo) - pad, int(hi) + pad
    if args.xlim:
        x0, x1 = args.xlim

    for p in (mark or ()):                       # SNP-mask positions as faint lines
        if x0 <= p + 1 <= x1:
            ax.axvline(p + 1, color="orange", lw=0.8, alpha=0.45, zorder=0)
    reads_by_umi = getattr(args, "_cell_reads", {}).get(gene) if args.reads_in_rows else None
    for i, u in enumerate(order, 1):
        ax.axhline(i, color="0.88", lw=0.5, ls=(0, (1, 2)), zorder=0)
        if reads_by_umi is not None and reads_by_umi.get(u):
            # one thin line per read inside the row of its transcript (UMI)
            rs = reads_by_umi[u]
            k = len(rs)
            band = 0.8
            for j, (r0, r1, uri, tc) in enumerate(rs):
                y = i - band / 2 + band * (j + 0.5) / k if k > 1 else i
                ax.plot([r0 + 1, r1], [y, y], ls=(0, (1, 1.5)), color="0.6", lw=0.6, zorder=0)
                ab = [p + 1 for p in uri if p not in tc or (p in snps and not args.show_masked)]
                ax.scatter(ab, [y] * len(ab), s=args.dot_size * 0.6, facecolors="0.86",
                           edgecolors="0.55", linewidths=0.3, zorder=1)
                xs = [p for p in tc if p not in snps]
                if xs:
                    ax.scatter([p + 1 for p in xs], [y] * len(xs),
                               c=[conv_n[u].get(p, 1) for p in xs], cmap=CMAP, norm=norm,
                               marker="x", s=args.dot_size * 1.8, linewidths=1.1, zorder=3)
                if args.show_masked:
                    m = [p + 1 for p in tc if p in snps]
                    if m:
                        ax.scatter(m, [y] * len(m), marker="s", s=args.dot_size * 1.4,
                                   facecolors="none", edgecolors="darkorange", linewidths=0.9,
                                   zorder=2)
            continue
        conv_sites = {p: c for p, c in conv_n[u].items()
                      if c >= args.min_conv_reads and (p not in snps or args.show_masked)}
        absent = [p + 1 for p in covered[u] if p not in conv_sites]
        ax.scatter(absent, [i] * len(absent), s=args.dot_size, facecolors="0.86",
                   edgecolors="0.55", linewidths=0.35, zorder=1)
        xs = [(p, c) for p, c in conv_sites.items() if p not in snps]
        if xs:
            ax.scatter([p + 1 for p, _ in xs], [i] * len(xs), c=[c for _, c in xs],
                       cmap=CMAP, norm=norm, marker="x", s=args.dot_size * 2.6,
                       linewidths=1.4, zorder=3)
        if args.show_masked:
            m = [p + 1 for p in conv_sites if p in snps]
            if m:
                ax.scatter(m, [i] * len(m), marker="s", s=args.dot_size * 2,
                           facecolors="none", edgecolors="darkorange", linewidths=1, zorder=2)

    ax.set_xlim(x0, x1)
    ax.set_ylim(n + 0.7, 0.3)
    if ORDER == "labeled_first" and order is order0:   # divider between the two groups
        nl = sum(1 for u in order if u in lab)
        if 0 < nl < n:
            ax.axhline(nl + 0.5, color="0.35", lw=0.9, ls="--", zorder=4)
        trans = matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData)
        if nl:
            ax.text(1.005, (nl + 1) / 2, f"labeled\n({nl})", transform=trans, ha="left",
                    va="center", fontsize=7, color="#1f1fff", clip_on=False)
        if n - nl:
            ax.text(1.005, (nl + 1 + n) / 2, f"unlabeled\n({n - nl})", transform=trans,
                    ha="left", va="center", fontsize=7, color="0.4", clip_on=False)
    ax.set_yticks([1, n] if n > 1 else [1])
    ax.tick_params(axis="y", length=3)
    n_all = len(lab) + len(unl)
    if n < n_all:   # --max-umis cut: labeled are shown first, so give the totals too
        ax.set_ylabel(f"{n} of {n_all} transcripts shown\n"
                      f"(all: {len(lab)} labeled, {len(unl)} unlabeled)", fontsize=9)
    else:
        ax.set_ylabel(f"{n} transcripts ({len(lab)} labeled, {len(unl)} unlabeled)", fontsize=9)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(6, integer=True))
    ax.grid(axis="x", color="0.85", lw=1.0, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel(f"{chrom} position", fontsize=9, labelpad=2)
    if title is None:
        title = f"$\\it{{{gene}}}$" + (f" ({label})" if label else "")
    ax.set_title(title, fontsize=11)
    for sp in ax.spines.values():
        sp.set_linewidth(0.8)
    return len(lab), len(unl)


def gtf_gene_names(gtf, contigs):
    names = set()
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as f:
        for line in f:
            if line.startswith("#") or "\tgene\t" not in line or 'gene_name "' not in line:
                continue
            c = line.split("\t", 1)[0]
            if c in contigs:
                names.add(line.split('gene_name "')[1].split('"')[0])
    return names


def choose_genes_by_ntr(a, contigs):
    """Pick high/low NTR genes (and a cell) from a dynast count adata.h5ad.
    NTR = labeled UMIs / all UMIs (layers labeled_TC/unlabeled_TC or X_l_TC/X_n_TC)."""
    if not a.adata:
        sys.exit("ERROR: --high-ntr-genes / --low-ntr-genes / --list-genes need --adata")
    import anndata
    import pandas as pd
    import scipy.sparse as sp
    ad = anndata.read_h5ad(a.adata)
    lab = unl = None
    for l1, l2 in (("labeled_TC", "unlabeled_TC"), ("X_l_TC", "X_n_TC")):
        if l1 in ad.layers and l2 in ad.layers:
            lab, unl = sp.csc_matrix(ad.layers[l1]), sp.csc_matrix(ad.layers[l2])
            break
    if lab is None:
        sys.exit(f"ERROR: no labeled/unlabeled TC layers in {a.adata}")
    names = (ad.var["gene_name"].astype(str).values if "gene_name" in ad.var
             else ad.var_names.astype(str).values)
    tot = lab + unl
    df = pd.DataFrame({
        "gene": names,
        "UMIs": np.asarray(tot.sum(0)).ravel().astype(int),
        "labeled": np.asarray(lab.sum(0)).ravel().astype(int),
        "cells": np.asarray((tot > 0).sum(0)).ravel().astype(int)})
    df["NTR"] = np.where(df.UMIs > 0, df.labeled / df.UMIs.clip(lower=1), 0.0)
    df["idx"] = np.arange(len(df))
    in_gtf = gtf_gene_names(a.gtf, contigs)
    keep = (df.UMIs >= a.min_gene_umis) & (df.cells >= a.min_cells) & df.gene.isin(in_gtf)
    if a.no_mt_ribo:
        keep &= ~df.gene.str.match(r"^(MT-|RPL|RPS)")
    keep &= ~df.gene.isin(a.genes)
    pool = df[keep].drop_duplicates("gene")
    pool = pool.assign(pctl=pool.NTR.rank(pct=True) * 100)     # NTR percentile among filtered
    # data-driven defaults (tuned on iPSC day 0: overall NTR ~29%, NTR falls with expression)
    if a.high_ntr_min is None and a.high_ntr_percentile is None:
        a.high_ntr_percentile = 90.0
    if a.low_ntr_max is None and a.low_ntr_percentile is None:
        a.low_ntr_percentile = 10.0
    if len(pool):
        q10, q90 = (pool.NTR.quantile([0.1, 0.9]) * 100).values
        if a.min_cell_ntr is None:
            a.min_cell_ntr = round(2 / 3 * q90, 1)
        if a.max_cell_ntr is None:
            a.max_cell_ntr = round(q10, 1)
    print(f"Settings: genes >= {a.min_gene_umis} UMIs, >= {a.min_cells} cells; high NTR: "
          f">= {a.min_gene_labeled} labeled"
          + (f", NTR >= {a.high_ntr_min:g}%" if a.high_ntr_min is not None else "")
          + (f", pctl >= {a.high_ntr_percentile:g}" if a.high_ntr_percentile else "")
          + "; low NTR:"
          + (f" NTR <= {a.low_ntr_max:g}%" if a.low_ntr_max is not None else "")
          + (f" pctl <= {a.low_ntr_percentile:g}" if a.low_ntr_percentile is not None else "")
          + f"; ranked by {a.rank_by}; cell: >= {a.min_cell_umis} UMIs per gene, "
          f">= {a.min_cell_labeled} labeled + NTR >= {a.min_cell_ntr}% (high), "
          f"NTR <= {a.max_cell_ntr}% (low)")
    hi_pool = pool[pool.labeled >= a.min_gene_labeled]
    lo_pool = pool
    if a.high_ntr_min is not None:
        hi_pool = hi_pool[hi_pool.NTR * 100 >= a.high_ntr_min]
    if a.high_ntr_percentile:
        hi_pool = hi_pool[hi_pool.pctl >= a.high_ntr_percentile]
    if a.low_ntr_max is not None:
        lo_pool = lo_pool[lo_pool.NTR * 100 <= a.low_ntr_max]
    if a.low_ntr_percentile is not None:
        lo_pool = lo_pool[lo_pool.pctl <= a.low_ntr_percentile]
    if a.rank_by == "umis":
        hi_pool = hi_pool.sort_values("UMIs", ascending=False)
        lo_pool = lo_pool.sort_values("UMIs", ascending=False)
    else:
        hi_pool = hi_pool.sort_values("NTR", ascending=False)
        lo_pool = lo_pool.sort_values("NTR")
    # thresholds given without a count: plot one gene of that kind
    if not a.high_ntr and not a.low_ntr and not a.genes and not a.list_genes:
        a.high_ntr = a.low_ntr = 1          # --adata alone: one high and one low NTR gene
    overall = df.labeled.sum() / max(df.UMIs.sum(), 1)
    print(f"NTR from {a.adata}: overall {overall:.2%}; {len(pool):,} genes pass filters "
          f"(>= {a.min_gene_umis} UMIs, >= {a.min_cells} cells, in GTF)")
    fmt = {"NTR": "{:.1%}".format, "pctl": "{:.0f}".format}
    cols = ["gene", "UMIs", "labeled", "cells", "NTR", "pctl"]
    q = pool.NTR.quantile([0.1, 0.5, 0.9]) * 100
    print(f"NTR of filtered genes: 10th pct {q.iloc[0]:.1f}%, median {q.iloc[1]:.1f}%, "
          f"90th pct {q.iloc[2]:.1f}%   (pctl column = NTR percentile)")
    show = a.list_genes or 10
    crit_hi = [f">= {a.min_gene_labeled} labeled UMIs"]
    if a.high_ntr_min is not None:
        crit_hi.append(f"NTR >= {a.high_ntr_min:g}%")
    if a.high_ntr_percentile:
        crit_hi.append(f"pctl >= {a.high_ntr_percentile:g}")
    crit_lo = []
    if a.low_ntr_max is not None:
        crit_lo.append(f"NTR <= {a.low_ntr_max:g}%")
    if a.low_ntr_percentile is not None:
        crit_lo.append(f"pctl <= {a.low_ntr_percentile:g}")
    print(f"\nHigh-NTR candidates ({', '.join(crit_hi)}; {len(hi_pool):,} genes, "
          f"ranked by {a.rank_by}):")
    print(hi_pool[cols].head(show).to_string(index=False, formatters=fmt))
    print(f"\nLow-NTR candidates ({', '.join(crit_lo) or 'no NTR limit'}; {len(lo_pool):,} genes, "
          f"ranked by {a.rank_by}):")
    print(lo_pool[cols].head(show).to_string(index=False, formatters=fmt))
    if a.list_genes:
        sys.exit(0)

    manual = list(a.genes)
    picked_hi = list(hi_pool.gene.head(a.high_ntr)) if a.high_ntr else []
    picked_lo = [g for g in lo_pool.gene if g not in picked_hi][: a.low_ntr] if a.low_ntr else []
    if a.high_ntr and len(picked_hi) < a.high_ntr:
        print(f"NOTE: only {len(picked_hi)} high-NTR gene(s) pass the filters; relax "
              "--high-ntr-percentile / --high-ntr-min / --min-gene-labeled / --min-gene-umis")
    if a.low_ntr and len(picked_lo) < a.low_ntr:
        print(f"NOTE: only {len(picked_lo)} low-NTR gene(s) pass the filters; relax "
              "--low-ntr-percentile / --low-ntr-max / --min-gene-umis")
    a.genes = picked_hi + manual + picked_lo
    if not a.labels:
        ntr = dict(zip(df.gene, df.NTR))
        a.labels = ([f"high NTR, {ntr[g]:.1%}" for g in picked_hi] +
                    [f"NTR {ntr[g]:.1%}" if g in ntr else "" for g in manual] +
                    [f"low NTR, {ntr[g]:.1%}" for g in picked_lo])
    print(f"\nGenes to plot: {', '.join(a.genes)}")

    # cell choice: every gene >= min_cell_umis UMIs, high-NTR genes >= min_cell_labeled labeled
    if a.cell or a.cell_mode != "same" or a.before_after or len(a.genes) == 0:
        return
    gi = {g: int(df.idx[df.gene == g].iloc[0]) for g in a.genes if (df.gene == g).any()}
    if len(gi) != len(a.genes):
        return
    T = np.column_stack([np.asarray(tot[:, gi[g]].todense()).ravel() for g in a.genes])
    L = np.column_stack([np.asarray(lab[:, gi[g]].todense()).ravel() for g in a.genes])
    ok = (T >= a.min_cell_umis).all(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(T > 0, L / np.maximum(T, 1), 0.0) * 100
    for j, g in enumerate(a.genes):
        if g in picked_hi:
            ok &= L[:, j] >= a.min_cell_labeled
            if a.min_cell_ntr is not None:
                ok &= R[:, j] >= a.min_cell_ntr
        if g in picked_lo and a.max_cell_ntr is not None:
            ok &= R[:, j] <= a.max_cell_ntr
    if not ok.any():
        print(f"NOTE: no cell passes the cell filters (>= {a.min_cell_umis} UMIs of every gene, "
              f">= {a.min_cell_labeled} labeled for high-NTR genes"
              + (f", cell NTR >= {a.min_cell_ntr:g}% (high)" if a.min_cell_ntr is not None else "")
              + (f", cell NTR <= {a.max_cell_ntr:g}% (low)" if a.max_cell_ntr is not None else "")
              + "); choosing the cell from the BAM instead")
        return
    score = np.where(ok, T.min(axis=1) + 0.001 * T.sum(axis=1), -1)
    best = np.argsort(-score)[: 10]
    print("\nBest cells (labeled / UMIs and NTR per gene):")
    for c in best:
        if score[c] < 0:
            break
        print("  " + str(ad.obs_names[c]) + "\t" + "  ".join(
            f"{g}={int(L[c, j])}/{int(T[c, j])} ({R[c, j]:.0f}%)" for j, g in enumerate(a.genes)))
    a.cell = str(ad.obs_names[best[0]])
    print(f"Using cell {a.cell}")


def rerun_self(drop, add):
    """Run this script again with the same options, minus `drop` (option names; their values
    are removed too) plus `add` (list of strings)."""
    import subprocess
    argv, out, i = sys.argv[1:], [], 0
    while i < len(argv):
        x = argv[i]
        name = x.split("=", 1)[0]
        if name in drop:
            i += 1
            if "=" not in x:                     # skip its value(s)
                while i < len(argv) and not argv[i].startswith("--"):
                    i += 1
            continue
        out.append(x)
        i += 1
    subprocess.run([sys.executable, os.path.abspath(__file__)] + out + add, check=False)


def choose_genes_by_mask(a, contigs):
    """--mask-genes N: genes whose labeled UMIs drop when the SNP mask is applied
    (dynast count without mask = --adata-unmasked, with the 0-based mask = --adata);
    filters --min-gene-umis, --min-removed, --min-removed-pct, --no-mt-ribo; order
    --mask-rank-by; --list-mask-genes N prints the table and stops."""
    if not (a.adata and a.adata_unmasked):
        sys.exit("ERROR: --mask-genes / --list-mask-genes need --adata (masked run) and "
                 "--adata-unmasked (run without mask)")
    import anndata
    import pandas as pd
    import scipy.sparse as sp

    def layers(path):
        ad = anndata.read_h5ad(path)
        for l1, l2 in (("labeled_TC", "unlabeled_TC"), ("X_l_TC", "X_n_TC")):
            if l1 in ad.layers and l2 in ad.layers:
                return ad, sp.csc_matrix(ad.layers[l1]), sp.csc_matrix(ad.layers[l2])
        sys.exit(f"ERROR: no labeled/unlabeled TC layers in {path}")
    ad1, lab1, unl1 = layers(a.adata)
    ad0, lab0, unl0 = layers(a.adata_unmasked)
    genes = ad0.var_names.intersection(ad1.var_names)
    i0, i1 = ad0.var_names.get_indexer(genes), ad1.var_names.get_indexer(genes)
    names = (ad0.var["gene_name"].astype(str).values if "gene_name" in ad0.var
             else ad0.var_names.astype(str).values)[i0]
    tot = np.asarray((lab0 + unl0).sum(0)).ravel()[i0]
    l0 = np.asarray(lab0.sum(0)).ravel()[i0]
    l1 = np.asarray(lab1.sum(0)).ravel()[i1]
    df = pd.DataFrame({"gene": names, "UMIs": tot.astype(int), "labeled_no_mask": l0.astype(int),
                       "labeled_mask": l1.astype(int)})
    df["removed"] = df.labeled_no_mask - df.labeled_mask
    df["NTR_no_mask"] = df.labeled_no_mask / df.UMIs.clip(lower=1)
    df["NTR_mask"] = df.labeled_mask / df.UMIs.clip(lower=1)
    df["removed_pct"] = np.where(df.labeled_no_mask > 0,
                                 100 * df.removed / df.labeled_no_mask.clip(lower=1), 0.0)
    in_gtf = gtf_gene_names(a.gtf, contigs)
    keep = ((df.UMIs >= a.min_gene_umis) & (df.removed >= max(a.min_removed, 1))
            & df.gene.isin(in_gtf))
    if a.min_removed_pct is not None:
        keep &= df.removed_pct >= a.min_removed_pct
    if a.no_mt_ribo:
        keep &= ~df.gene.str.match(r"^(MT-|RPL|RPS)")
    keep &= ~df.gene.isin(a.genes)
    key = {"removed": ["removed", "UMIs"], "pct": ["removed_pct", "removed"],
           "umis": ["UMIs", "removed"], "ntr": ["NTR_mask", "removed"]}[a.mask_rank_by]
    pool = df[keep].sort_values(key, ascending=False).drop_duplicates("gene")
    print(f"Settings: genes >= {a.min_gene_umis} UMIs, >= {a.min_removed} labeled UMIs removed"
          + (f", >= {a.min_removed_pct:g}% of labeled removed" if a.min_removed_pct is not None
             else "") + f"; ranked by {a.mask_rank_by}")
    print(f"Genes changed by the mask: {len(pool):,}; labeled UMIs removed in total: "
          f"{int(df.removed.clip(lower=0).sum()):,} of {int(df.labeled_no_mask.sum()):,}")
    fmt = {"NTR_no_mask": "{:.1%}".format, "NTR_mask": "{:.1%}".format,
           "removed_pct": "{:.0f}%".format}
    show = a.list_mask_genes or max(a.mask_genes or 0, 10)
    print(pool.head(show)[["gene", "UMIs", "labeled_no_mask", "labeled_mask", "removed",
                           "removed_pct", "NTR_no_mask", "NTR_mask"]]
          .to_string(index=False, formatters=fmt))
    if a.list_mask_genes:
        sys.exit(0)
    if pool.empty:
        sys.exit("No gene is changed by the mask with these filters "
                 "(lower --min-gene-umis / --min-removed / --min-removed-pct)")
    picked = list(pool.gene.head(a.mask_genes))
    if len(picked) < a.mask_genes:
        print(f"NOTE: only {len(picked)} gene(s) pass the filters")
    a.genes = list(a.genes) + picked
    print(f"Genes to plot: {', '.join(a.genes)}")

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g0 = ap.add_argument_group("Input (required)")
    g0.add_argument("--bam", required=True, help="PRE-consensus BAM (indexed)")
    g0.add_argument("--fasta", required=True)
    g0.add_argument("--gtf", required=True)
    g0.add_argument("--genes", nargs="+", default=[],
                    help="gene names, top to bottom (or pick them with --high-ntr-genes/--low-ntr-genes)")
    gn = ap.add_argument_group(
        "Gene choice by NTR (new-to-total ratio = labeled / all UMIs, from a dynast count "
        "adata.h5ad; use the masked run)")
    gn.add_argument("--adata", help="dynast count adata.h5ad (e.g. 02_dynast_count/TC_rm_bg_0based)")
    gn.add_argument("--high-ntr-genes", "--high-ntr", dest="high_ntr", type=int, default=0,
                    metavar="N",
                    help="number of high-NTR genes to plot (fast turnover; top panels). "
                         "Old name: --high-ntr")
    gn.add_argument("--low-ntr-genes", "--low-ntr", dest="low_ntr", type=int, default=0,
                    metavar="N",
                    help="number of low-NTR genes to plot (stable; bottom panels). "
                         "Old name: --low-ntr")
    gn.add_argument("--min-gene-umis", type=int, default=500,
                    help="gene filter: >= this many UMIs over all cells")
    gn.add_argument("--min-gene-labeled", type=int, default=50,
                    help="high-NTR filter: >= this many labeled UMIs over all cells")
    gn.add_argument("--min-cells", type=int, default=20,
                    help="gene filter: expressed (>= 1 UMI) in >= this many cells")
    gn.add_argument("--no-mt-ribo", action="store_true", help="skip MT-, RPL*, RPS* genes")
    gn.add_argument("--min-cell-umis", type=int, default=5,
                    help="cell choice: every plotted gene needs >= this many UMIs in the cell")
    gn.add_argument("--min-cell-labeled", type=int, default=2,
                    help="cell choice: every high-NTR gene needs >= this many labeled UMIs "
                         "in the cell")
    gn.add_argument("--high-ntr-min", type=float, metavar="PCT",
                    help="high-NTR genes need NTR >= PCT %% (e.g. 20)")
    gn.add_argument("--low-ntr-max", type=float, metavar="PCT",
                    help="low-NTR genes need NTR <= PCT %% (e.g. 2)")
    gn.add_argument("--high-ntr-percentile", type=float, metavar="P",
                    help="high-NTR genes need NTR >= the P-th percentile of the filtered genes "
                         "(default 90 = top 10%% when no --high-ntr-min; 0 = off)")
    gn.add_argument("--low-ntr-percentile", type=float, metavar="P",
                    help="low-NTR genes need NTR <= the P-th percentile "
                         "(default 10 = bottom 10%% when no --low-ntr-max; 100 = off)")
    gn.add_argument("--rank-by", choices=["ntr", "umis"], default="umis",
                    help="among genes passing the filters, pick by NTR (most extreme first) or "
                         "by UMIs (best expressed first)")
    gn.add_argument("--min-cell-ntr", type=float, metavar="PCT",
                    help="cell choice: every high-NTR gene needs NTR >= PCT %% in the cell "
                         "(default: 2/3 of the gene-level 90th-percentile NTR; 0 = off)")
    gn.add_argument("--max-cell-ntr", type=float, metavar="PCT",
                    help="cell choice: every low-NTR gene needs NTR <= PCT %% in the cell "
                         "(default: the gene-level 10th-percentile NTR; 100 = off)")
    gn.add_argument("--adata-unmasked",
                    help="dynast count adata.h5ad WITHOUT the SNP mask (e.g. 02_dynast_count/TC); "
                         "with --mask-genes")
    gn.add_argument("--mask-genes", type=int, metavar="N",
                    help="before-after plots: pick the N genes whose labeled UMIs drop most with "
                         "the mask (needs --adata = masked run and --adata-unmasked)")
    gn.add_argument("--list-mask-genes", type=int, metavar="N",
                    help="print the N genes the mask changes most (after filters) and exit")
    gn.add_argument("--min-removed", type=int, default=5,
                    help="--mask-genes filter: >= this many labeled UMIs removed by the mask")
    gn.add_argument("--min-removed-pct", type=float, metavar="PCT",
                    help="--mask-genes filter: >= PCT %% of the gene's labeled UMIs removed")
    gn.add_argument("--mask-rank-by", choices=["removed", "pct", "umis", "ntr"], default="removed",
                    help="--mask-genes order: most labeled UMIs removed, highest %% removed, "
                         "best expressed, or highest NTR after masking")
    gn.add_argument("--list-genes", type=int, metavar="N",
                    help="print the N highest- and lowest-NTR genes (after filters) and exit")
    g1 = ap.add_argument_group("Plot type (choose one; default = standard transcript plot)")
    g1.add_argument("--before-after", action="store_true",
                    help="ONE gene, two panels: counted without the SNP mask (top) and with it "
                         "(bottom), same cell and same transcript order")
    g1.add_argument("--scnt-layout", "--paper", dest="paper", action="store_true",
                    help="scNT-seq figure layout (Qiu et al. 2020, ACTG1 3'UTR figure) for ONE "
                         "gene: gene bar with strand arrows, one row per transcript, coverage "
                         "colours 1/2/3/>=4, a highlighted transcript and an inset with all its "
                         "reads (use the pre-consensus BAM). --paper is the old name")
    g1.add_argument("--compare-consensus",
                    help="ONE gene, 2x2 panels from the same cell and transcripts: --bam "
                         "(pre-consensus) vs this consensus BAM, each without and with the SNP "
                         "mask (needs --snp-csv). With --before-after, the site/cell/window "
                         "are chosen as usual")
    g1.add_argument("--umi-compare", action="store_true",
                    help="ONE UMI in two panels: its reads in --bam (pre-consensus) and its "
                         "record(s) in --compare-consensus (consensus.bam). UMI from "
                         "--umi-reads, or the labeled UMI with the most reads")
    g1.add_argument("--umi-reads", nargs="?", const="auto",
                    help="add a panel with every read of ONE transcript (UMI) on its own line, "
                         "like the scNT-seq inset (use the pre-consensus BAM). Give a UMI, or no "
                         "value for the UMI with the most reads (labeled UMIs preferred)")
    g1.add_argument("--umi-only", action="store_true",
                    help="with --umi-reads: draw ONLY the reads of that UMI (no transcript panel)")
    g1.add_argument("--list-cells", type=int, metavar="N",
                    help="only print the N cells with the most UMIs for each gene (one line per "
                         "cell: gene, cell, UMIs) and exit; use it to loop over top cells")
    g2 = ap.add_argument_group("Cell and transcript selection")
    g2.add_argument("--cell", help="cell barcode; default: best cell across all genes")
    g2.add_argument("--top-cells", type=int, metavar="N",
                    help="make one figure for each of the N best cells (most UMIs in ALL listed "
                         "genes); files get _top1, _top2, ... and the barcode")
    g2.add_argument("--cell-mode", choices=["same", "each"], default="same",
                    help="same: one cell for all panels (as in the paper); "
                         "each: the cell with most UMIs for each gene separately")
    g2.add_argument("--highlight-umis", nargs="+",
                    help="--compare-consensus: UMI(s) to highlight in every panel "
                         "(default: the labeled transcript with the most reads)")
    g2.add_argument("--no-highlight", action="store_true",
                    help="--compare-consensus: do not highlight any transcript")
    g2.add_argument("--max-umis", type=int, default=60, help="max transcripts per panel")
    g2.add_argument("--subsample", choices=["labeled_first", "proportional", "first"],
                    default="labeled_first",
                    help="if a cell has more than --max-umis transcripts: keep all labeled and a "
                         "random set of unlabeled (labeled_first), a random set keeping the "
                         "labeled fraction (proportional), or the leftmost ones (first)")
    g2.add_argument("--group-by-label", choices=["yes", "no"],
                    help="yes: labeled transcripts on top, unlabeled below, with a divider line and "
                         "group names; no: all transcripts by position (one diagonal). "
                         "Short for --order labeled_first / position (default: yes)")
    g2.add_argument("--order", choices=["labeled_first", "position"], default="labeled_first",
                    help="row order: labeled transcripts on top, each group by position "
                         "(labeled_first, as in the scNT-seq figure), or all transcripts by "
                         "position (position)")
    g3 = ap.add_argument_group("SNP mask and negative controls")
    g3.add_argument("--snp-csv", help="SNP mask (corrected, 0-based)")
    g3.add_argument("--snp-coord", choices=["0based", "1based"], default="0based")
    g3.add_argument("--show-masked", action="store_true",
                    help="draw T>C at masked SNPs as orange squares instead of grey")
    g3.add_argument("--only-masked", action="store_true",
                    help="--before-after: draw only the panel counted with the corrected mask "
                         "(site, cell and window are still chosen automatically)")
    g3.add_argument("--control-shift", type=int, nargs="*", default=[],
                    help="--before-after negative control: add a panel counted with the SNP mask "
                         "shifted by N bp (e.g. 1 = a mask built from VCF POS without -1, as "
                         "dynast would read it). Several values allowed, e.g. 1 -1")
    g3.add_argument("--control-mask",
                    help="--before-after negative control: add a panel counted with another mask "
                         "file, read as 0-based (e.g. an unrelated donor's mask)")
    g3.add_argument("--control-label", default="control mask",
                    help="panel title for --control-mask")
    g3.add_argument("--old-snp-csv", help=argparse.SUPPRESS)
    g4 = ap.add_argument_group("Pre-consensus vs consensus comparison (with --compare-consensus)")
    g4.add_argument("--compare-layout", choices=["inset", "side", "mask"], default="inset",
                    help="--scnt-layout with --compare-consensus: 'inset' adds the consensus record "
                         "below the reads in the inset; 'side' draws the whole paper figure "
                         "twice, side by side (pre-consensus BAM | consensus.bam), same cell, "
                         "rows and highlighted UMI")
    g4.add_argument("--compare-panels", choices=["2x2", "masked", "nomask"], default="2x2",
                    help="--compare-consensus: 2x2 (both BAMs x without/with mask), or two "
                         "panels side by side, pre-consensus vs consensus, with the SNP mask "
                         "(masked) or without (nomask)")
    g4.add_argument("--consensus-min-bq", type=int, default=28,
                    help="--compare-consensus: base-quality cutoff for the consensus BAM "
                         "(dynast counts quality > 27); --min-bq applies to --bam")
    g5 = ap.add_argument_group("Genomic window")
    g5.add_argument("--xlim", type=int, nargs=2, help="1-based window (zoom on a SNP)")
    g5.add_argument("--utr3", action="store_true",
                    help="show the 3'UTR of the gene's canonical transcript (one gene only), "
                         "like the scNT-seq ACTG1 3'UTR figure")
    g5.add_argument("--site", type=int,
                    help="--before-after: 1-based position of the masked SNP to show "
                         "(default: the site where the mask removes the most T>C)")
    g5.add_argument("--window", type=int, default=150,
                    help="--before-after: bp shown on each side of the site (default 150)")
    g5.add_argument("--full-range", action="store_true",
                    help="default window: show every covered position (by default the window "
                         "covers the central part of the coverage, see --range-pct)")
    g5.add_argument("--range-pct", type=float, default=2.5,
                    help="default window: drop this percent of covered positions on each side "
                         "(default 2.5)")
    g6 = ap.add_argument_group("Appearance and output")
    g6.add_argument("--labels", nargs="*", default=[], help="subtitle per gene (optional)")
    g6.add_argument("--sample-label", default="", help="--scnt-layout: panel title, e.g. 'Day 0'")
    g6.add_argument("--label-umis", action="store_true",
                    help="--scnt-layout / --compare-consensus: write each row's UMI sequence on the "
                         "y-axis; the highlighted UMI in bold red")
    g6.add_argument("--umi-fontsize", type=float, default=6)
    g6.add_argument("--reads-in-rows", action="store_true",
                    help="each row is one transcript (UMI) and shows each of its reads as a thin "
                         "line inside the row (use the pre-consensus BAM)")
    g6.add_argument("--row-height", type=float,
                    help="inches per transcript row (default 0.11, or 0.22 with --reads-in-rows)")
    g6.add_argument("--vmax", type=float, default=16, help="colour scale max (paper: 16)")
    g6.add_argument("--dot-size", type=float, default=7)
    g6.add_argument("--width", type=float, default=7.5)
    g6.add_argument("--out",
                    help="output file; default: <plot name>_group-<yes|no>_mask-<yes|no>_<genes>.pdf "
                         "in --outdir, e.g. basic_group-yes_mask-yes_SOX2-ACTB.pdf (the cell "
                         "barcode is always appended)")
    g6.add_argument("--outdir", default=".", help="folder for the automatic file name")
    g6.add_argument("--format", default="pdf", choices=["pdf", "png", "svg"],
                    help="file type for the automatic file name")
    g7 = ap.add_argument_group("Read filters")
    g7.add_argument("--strand", choices=["forward", "reverse", "unstranded"], default="forward")
    g7.add_argument("--cb-tag", default="CB")
    g7.add_argument("--ub-tag", default="UB")
    g7.add_argument("--min-bq", type=int, default=20)
    g7.add_argument("--min-mapq", type=int, default=255)
    g7.add_argument("--min-conv-reads", type=int, default=1)
    a = ap.parse_args()
    global SUBSAMPLE, ORDER
    if a.group_by_label:
        a.order = "labeled_first" if a.group_by_label == "yes" else "position"
    SUBSAMPLE, ORDER = a.subsample, a.order

    bam = pysam.AlignmentFile(a.bam, "rb")
    fa = pysam.FastaFile(a.fasta)
    if a.mask_genes or a.list_mask_genes:
        a.before_after = True
        choose_genes_by_mask(a, set(bam.references))
    elif a.adata or a.high_ntr or a.low_ntr or a.list_genes:
        choose_genes_by_ntr(a, set(bam.references))
    if (a.before_after or (a.paper and a.compare_layout == "mask")) and len(a.genes) > 1:
        # one mask figure per gene
        for g in a.genes:
            print(f"\n===== {g}")
            extra = []
            if a.out:                                # keep user names apart per gene
                root, ext = os.path.splitext(a.out)
                extra = ["--out", f"{root}_{g}{ext}"]
            rerun_self({"--genes", "--mask-genes", "--list-mask-genes", "--adata",
                        "--adata-unmasked", "--labels", "--out"},
                       ["--genes", g, "--before-after"] + extra)
        return
    if not a.genes:
        sys.exit("ERROR: give --genes, or --adata with --high-ntr-genes/--low-ntr-genes")
    labels = a.labels + [""] * (len(a.genes) - len(a.labels))
    locs = find_genes(a.gtf, a.genes, set(bam.references))
    global PLOT_NAME
    PLOT_NAME = plot_name(a)
    if not a.out:
        os.makedirs(a.outdir, exist_ok=True)
        a.out = os.path.join(a.outdir, f"{PLOT_NAME}_{'_'.join(name_tags(a))}_"
                                       f"{'-'.join(a.genes)}.{a.format}")
    else:
        a.out = add_tags(a.out, name_tags(a))     # your --out name + group-/mask- tags
    if os.path.dirname(a.out):
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
    print(f"Plot type: {PLOT_NAME}")
    for g in a.genes:
        c, s, e, st = locs[g]
        print(f"{g}: {c}:{s + 1:,}-{e:,} ({st})")
    if a.xlim:
        c, s, e, _ = locs[a.genes[0]]
        if a.xlim[1] < s + 1 or a.xlim[0] > e:
            sys.exit(f"ERROR: --xlim {a.xlim[0]:,}-{a.xlim[1]:,} is outside {a.genes[0]} "
                     f"({c}:{s + 1:,}-{e:,}); use --utr3 or coordinates inside the gene")
    if a.utr3:
        if len(a.genes) != 1:
            sys.exit("ERROR: --utr3 works with one gene")
        u0, u1 = utr3_span(a.gtf, a.genes[0])
        a.xlim = [u0, u1]
        print(f"{a.genes[0]} 3'UTR: {locs[a.genes[0]][0]}:{u0:,}-{u1:,} ({u1 - u0 + 1} bp)")

    if a.top_cells and not a.cell and not a.before_after and \
            not (a.paper and a.compare_layout == "mask"):
        import subprocess
        per_gene = {g: umis_per_cell(bam, locs[g], a) for g in a.genes}
        common = set.intersection(*[set(per_gene[g]) for g in a.genes])
        ranked = sorted(common, key=lambda c: (min(per_gene[g][c] for g in a.genes),
                                               sum(per_gene[g][c] for g in a.genes)),
                        reverse=True)[: a.top_cells]
        if not ranked:
            sys.exit("ERROR: no cell has UMIs in all listed genes")
        argv = sys.argv[1:]
        if "--top-cells" in argv:                 # drop "--top-cells N"
            i = argv.index("--top-cells"); del argv[i:i + 2]
        argv = [x for x in argv if not x.startswith("--top-cells=")]
        if "--out" in argv:                       # set --out per rank
            i = argv.index("--out"); del argv[i:i + 2]
        root, ext = os.path.splitext(a.out)
        print(f"Top {len(ranked)} cells (UMIs per gene):")
        for r, c in enumerate(ranked, 1):
            print(f"  #{r} {c}  " + "  ".join(f"{g}={per_gene[g][c]}" for g in a.genes))
        for r, c in enumerate(ranked, 1):
            print(f"\n===== #{r} cell {c}")
            subprocess.run([sys.executable, os.path.abspath(__file__)] + argv +
                           ["--cell", c, "--out", f"{root}_top{r}{ext}"], check=False)
        return

    if a.list_cells:
        for g in a.genes:
            d = umis_per_cell(bam, locs[g], a)
            for cb, n in sorted(d.items(), key=lambda kv: -kv[1])[: a.list_cells]:
                print(f"TOPCELL\t{g}\t{cb}\t{n}")
        return

    per_gene = None
    inset_mask = a.paper and a.compare_layout == "mask"
    if a.before_after or inset_mask:
        # mask figures (transcript panels or inset layout): zoom on the masked SNP that
        # removes most T>C, in the cell with the most molecules carrying it
        if len(a.genes) != 1 or not a.snp_csv:
            sys.exit("ERROR: mask comparison needs exactly one gene and --snp-csv")
        g = a.genes[0]
        c, s0, e0, st = locs[g]
        gsnps = {p for p in load_snps(a.snp_csv, c, a.snp_coord) if s0 <= p < e0}
        print(f"{len(gsnps)} masked SNP sites in {g}; scanning reads for T>C at these sites ...")
        sites = scan_removed(bam, fa, locs[g], gsnps, a)
        if not sites:
            print(f"NOTE: no T>C at any masked SNP site in {g} for this sample -- the mask does "
                  f"not change {g} here; plotting anyway (panels will be identical)")
            if not a.cell:
                d = umis_per_cell(bam, locs[g], a)
                if not d:
                    sys.exit(f"ERROR: no reads with cell barcode and UMI for {g}")
                a.cell = max(d, key=d.get)
                print(f"Using the cell with the most {g} UMIs: {a.cell} ({d[a.cell]} UMIs)")
        if sites:
            ranked = sorted(sites.items(), key=lambda kv: -kv[1]["umis"])
            print("Masked sites with removed T>C (all cells):")
            print("  pos_1based\tUMIs_with_TC\treads_with_TC\tcells")
            for p, d in ranked[:10]:
                print(f"  {c}:{p + 1:,}\t{d['umis']}\t{d['reads']}\t{len(d['cells'])}")
            if a.site:
                site = a.site - 1
                if site not in sites:
                    sys.exit(f"ERROR: {c}:{a.site:,} has no removed T>C in this sample")
            elif a.cell:
                # given cell: use the site where THIS cell has the most variant molecules
                cb = norm_bc(a.cell)
                in_cell = [(p, d["cells"][cb][0]) for p, d in sites.items() if cb in d["cells"]]
                site = max(in_cell, key=lambda x: x[1])[0] if in_cell else None
                if site is None:
                    print(f"NOTE: cell {cb} has no T>C at any masked site of {g}; plotting its "
                          "covered region (panels will be identical)")
            else:
                site = ranked[0][0]
            if site is not None:
                best = sorted(sites[site]["cells"].items(), key=lambda kv: (-kv[1][0], -kv[1][1]))
                print(f"Cells at {c}:{site + 1:,} (UMIs with T>C / UMIs covering the site):")
                for cb, (k, n) in best[:10]:
                    print(f"  {cb}\t{k}/{n}")
                if a.top_cells and not a.cell:   # one figure per cell, same site
                    root, ext = os.path.splitext(a.out)
                    for r, (cb, (k, n)) in enumerate(best[: a.top_cells], 1):
                        print(f"\n===== #{r} cell {cb} ({k}/{n} molecules with T>C at the site)")
                        rerun_self({"--top-cells", "--out", "--cell", "--site", "--genes",
                                    "--mask-genes", "--adata", "--adata-unmasked"},
                                   ["--genes", g, "--before-after", "--cell", cb,
                                    "--site", str(site + 1), "--out", f"{root}_top{r}{ext}"])
                    return
                if not a.cell:
                    a.cell = best[0][0]
                if not a.xlim:
                    a.xlim = [site + 1 - a.window, site + 1 + a.window]
                print(f"Showing {c}:{site + 1:,} in cell {norm_bc(a.cell)}, window {a.xlim[0]:,}-{a.xlim[1]:,}")
    if a.cell:
        cells = {g: norm_bc(a.cell) for g in a.genes}
    else:
        per_gene = {g: umis_per_cell(bam, locs[g], a) for g in a.genes}
        if a.cell_mode == "same":
            c = pick_cell(per_gene, a.genes)
            cells = {g: c for g in a.genes}
        else:
            report_genes(per_gene, a.genes)
            for g in a.genes:
                if not per_gene[g]:
                    sys.exit(f"ERROR: no cell has UMIs for {g} in this BAM")
            cells = {g: max(per_gene[g], key=per_gene[g].get) for g in a.genes}
    for g in a.genes:
        print(f"{g}: cell {cells[g]}")

    data = {}
    for g in a.genes:
        c, s, e, st = locs[g]
        snps = set()
        if a.snp_csv:
            snps = {p for p in load_snps(a.snp_csv, c, a.snp_coord) if s <= p < e}
        cov, conv, _ = collect(bam, fa, c, s, e, st, cells[g], a)
        data[g] = (cov, conv, snps)

    if a.reads_in_rows:
        a._cell_reads = {g: collect_cell_reads(bam, fa, locs[g], cells[g], a) for g in a.genes}
        mx = max((len(v) for d in a._cell_reads.values() for v in d.values()), default=1)
        print(f"--reads-in-rows: up to {mx} reads per transcript"
              + ("  (1 per transcript: consensus BAM? use the pre-consensus BAM)" if mx <= 1 else ""))
    if a.xlim:   # keep only uridines inside the window; UMIs outside it drop out
        lo, hi = a.xlim[0] - 1, a.xlim[1]
        for g, (cov, conv, snps) in data.items():
            cov2 = defaultdict(set, {u: {p for p in ps if lo <= p < hi} for u, ps in cov.items()})
            conv2 = defaultdict(dict, {u: {p: c for p, c in cs.items() if lo <= p < hi}
                                       for u, cs in conv.items()})
            for u in list(conv2):
                conv2[u] = defaultdict(int, conv2[u])
            data[g] = (cov2, conv2, snps)

    if a.umi_compare:
        if len(a.genes) != 1 or not a.compare_consensus:
            sys.exit("ERROR: --umi-compare needs one gene and --compare-consensus <consensus.bam>")
        g = a.genes[0]
        cell = cells[g]
        snps = data[g][2]
        pre = collect_cell_reads(bam, fa, locs[g], cell, a)
        bq = a.min_bq
        a.min_bq = a.consensus_min_bq
        con = collect_cell_reads(pysam.AlignmentFile(a.compare_consensus, "rb"), fa, locs[g], cell, a)
        a.min_bq = bq
        if a.umi_reads and a.umi_reads != "auto":
            umi = a.umi_reads
        else:   # labeled UMI (T>C outside the mask) with the most reads
            umi = max(pre, key=lambda u: (any(p not in snps for r in pre[u] for p in r[3]),
                                          len(pre[u])), default=None)
        if not umi or umi not in pre:
            sys.exit(f"ERROR: UMI {umi} has no reads for {g} in cell {cell}")
        rp, rc = pre[umi], con.get(umi, [])
        print(f"UMI {umi} in cell {cell}: {len(rp)} pre-consensus read(s), "
              f"{len(rc)} consensus record(s)")
        npos = defaultdict(int)          # reads with T>C per position (pre-consensus)
        for r in rp:
            for p in r[3]:
                npos[p] += 1
        allr = rp + rc
        lo_ = min(r[0] for r in allr) + 1
        hi_ = max(r[1] for r in allr)
        pad = max(8, int(0.06 * (hi_ - lo_)))
        if a.xlim:
            lo_, hi_, pad = a.xlim[0], a.xlim[1], 0

        def draw(ax, reads, cov_of, title, line_col):
            for j, (r0, r1, uri, tc) in enumerate(reads, 1):
                ax.plot([r0 + 1, r1], [j, j], ls=(0, (1.5, 1.5)), color=line_col, lw=0.9, zorder=0)
                ab = [p + 1 for p in uri if p not in tc]
                ax.scatter(ab, [j] * len(ab), s=a.dot_size * 1.3, facecolors="0.88",
                           edgecolors="0.45", linewidths=0.4, zorder=1)
                xs = [p for p in tc if p not in snps]
                if xs:
                    ax.scatter([p + 1 for p in xs], [j] * len(xs), marker="x",
                               color=[PAPER_COLORS[min(cov_of(p), 4) - 1] for p in xs],
                               s=a.dot_size * 3.2, linewidths=1.7, zorder=3)
                ms = [p + 1 for p in tc if p in snps]
                if ms:
                    ax.scatter(ms, [j] * len(ms), marker="s", s=a.dot_size * 2.4,
                               facecolors="none", edgecolors="darkorange", linewidths=1.0,
                               zorder=2)
            k = max(len(reads), 1)
            ax.set_ylim(k + 0.7, 0.3)
            ax.set_yticks(range(1, len(reads) + 1))
            ax.set_xlim(lo_ - pad, hi_ + pad)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5, integer=True))
            ax.grid(axis="x", color="0.85", lw=1.0, zorder=0)
            ax.set_axisbelow(True)
            ax.set_title(title, fontsize=10)
            for sp in ax.spines.values():
                sp.set_linewidth(1.2)
            if not reads:
                ax.text(0.5, 0.5, "UMI not in this BAM (window / filters)", ha="center",
                        transform=ax.transAxes, color="0.4")

        hp = max(1.4, 0.3 * len(rp) + 0.9)
        hc = max(1.1, 0.3 * max(len(rc), 1) + 0.9)
        fig = plt.figure(figsize=(a.width, 1.6 + hp + hc + 0.9))
        gs = fig.add_gridspec(3, 1, height_ratios=[1.4, hp, hc],
                              hspace=0.9 / np.mean([1.4, hp, hc]))
        lax = fig.add_subplot(gs[0]); lax.axis("off")
        lax.text(0.5, 1.0, f"$\\it{{{g}}}$  cell {cell}  UMI {umi}", ha="center", va="top",
                 fontsize=11, color="#e8302a", family="monospace", transform=lax.transAxes)
        h1 = [Line2D([], [], marker="o", ls="", mfc="0.88", mec="0.45", ms=8, label="Absent"),
              Line2D([], [], marker="x", ls="", color=PAPER_COLORS[0], mew=2, ms=8, label="Present")]
        if snps:
            h1.append(Line2D([], [], marker="s", ls="", mfc="none", mec="darkorange", ms=8,
                             label="Masked SNP (not counted)"))
        l1 = lax.legend(handles=h1, title="T-to-C substitution at a uridine site:",
                        loc="upper center", bbox_to_anchor=(0.5, 0.82), frameon=False, fontsize=9,
                        title_fontsize=9, ncol=len(h1))
        lax.add_artist(l1)
        lax.legend(handles=[Line2D([], [], marker="s", ls="", color=c, ms=9,
                                   label=("≥4" if i == 4 else str(i)))
                            for i, c in enumerate(PAPER_COLORS, 1)],
                   title="Reads with T-to-C at the site:", loc="upper center",
                   bbox_to_anchor=(0.5, 0.35), frameon=False, fontsize=9, title_fontsize=9,
                   ncol=4, handletextpad=0.2, columnspacing=0.8)
        ax1 = fig.add_subplot(gs[1])
        draw(ax1, rp, lambda p: npos[p],
             f"Pre-consensus BAM: {len(rp)} read(s) with this UMI (base quality >= {a.min_bq})",
             "0.45")
        ax1.set_yticklabels([f"read {j}" for j in range(1, len(rp) + 1)], fontsize=8)
        ax2 = fig.add_subplot(gs[2], sharex=ax1)
        draw(ax2, rc, lambda p: 1,
             f"consensus.bam: {len(rc)} record(s) for this UMI "
             f"(base quality >= {a.consensus_min_bq})", "#d97706")
        ax2.set_yticklabels(["consensus"] * len(rc), fontsize=8, color="#d97706")
        ax2.set_xlabel(f"{locs[g][0]} position", fontsize=9)
        lab_pre = any(p not in snps for r in rp for p in r[3])
        lab_con = any(p not in snps for r in rc for p in r[3])
        print(f"labeled: pre-consensus {'yes' if lab_pre else 'no'}, "
              f"consensus {'yes' if lab_con else 'no'}")
        for p in sorted(set(npos) | {q for r in rc for q in r[3]}):
            inc = any(p in r[3] for r in rc)
            print(f"  {locs[g][0]}:{p + 1:,}  T>C in {npos.get(p, 0)}/{len(rp)} reads  "
                  f"consensus: {'T>C' if inc else '-'}" + ("  (masked SNP)" if p in snps else ""))
        label_cells(fig, [cell])
        a.out = out_with_cell(a.out, [cell])
        fig.savefig(a.out, dpi=300, bbox_inches="tight")
        print(f"Saved: {a.out}")
        return

    if a.compare_consensus and not a.paper:
        if len(a.genes) != 1:
            sys.exit("ERROR: --compare-consensus needs exactly one gene")
        if a.compare_panels != "nomask" and not a.snp_csv:
            sys.exit("ERROR: --compare-consensus with mask panels needs --snp-csv "
                     "(or use --compare-panels nomask)")
        g = a.genes[0]
        c_, s_, e_, st_ = locs[g]
        cov_p, conv_p, snps = data[g]
        bam_c = pysam.AlignmentFile(a.compare_consensus, "rb")
        bq_pre = a.min_bq
        a.min_bq = a.consensus_min_bq
        cov_c, conv_c, nr_c = collect(bam_c, fa, c_, s_, e_, st_, cells[g], a)
        a.min_bq = bq_pre
        if a.xlim:
            lo, hi = a.xlim[0] - 1, a.xlim[1]
            cov_c = defaultdict(set, {u: {p for p in ps if lo <= p < hi} for u, ps in cov_c.items()})
            conv_c = defaultdict(lambda: defaultdict(int),
                                 {u: defaultdict(int, {p: k for p, k in cs.items() if lo <= p < hi})
                                  for u, cs in conv_c.items()})
        # same transcripts and order in all panels: union of UMIs, ordered from the
        # pre-consensus BAM without mask (UMIs only in one BAM are appended)
        cov_all = defaultdict(set)
        for u, ps in list(cov_p.items()) + list(cov_c.items()):
            cov_all[u] |= ps
        order, _, _ = order_umis(cov_all, conv_p, set(), a.min_conv_reads, a.max_umis)
        # transcripts to highlight (red box + shading in every panel):
        #   --highlight-umis U1 U2 ...  or  --umi-reads U  -> those UMIs
        #   otherwise (default)                            -> the labeled transcript with the
        #                                                     most reads; --no-highlight: none
        reps = []
        if a.highlight_umis:
            reps = list(a.highlight_umis)
        elif a.umi_reads and a.umi_reads != "auto":
            reps = [a.umi_reads]
        elif not a.no_highlight:
            cr = collect_cell_reads(bam, fa, locs[g], cells[g], a)
            lab_p = {u for u in order if any(p not in snps for p in conv_p[u])}
            cand = [u for u in order if cr.get(u)]
            if cand:
                reps = [max(cand, key=lambda u: (u in lab_p, len(cr[u])))]
        missing = [u for u in reps if u not in order]
        if missing:
            print(f"NOTE: not among the shown transcripts (raise --max-umis?): {', '.join(missing)}")
        reps = [u for u in reps if u in order]
        for u in reps:
            print(f"Highlighted transcript: UMI {u} (row {order.index(u) + 1})")
        rep = reps[0] if reps else None
        norm = Normalize(vmin=1, vmax=a.vmax)
        nrow = len(order)
        h = max(2.0, (a.row_height or (0.22 if a.reads_in_rows else 0.11)) * nrow + 1.0)
        two = a.compare_panels != "2x2"
        fig = plt.figure(figsize=(a.width * 1.8, (h if two else 2 * h) + 1.0))
        if two:
            gs = fig.add_gridspec(2, 2, height_ratios=[0.6, h], hspace=1.0 / np.mean([0.6, h]),
                                  wspace=0.18)
        else:
            gs = fig.add_gridspec(3, 2, height_ratios=[0.6, h, h], hspace=0.5, wspace=0.18)
        lax = fig.add_subplot(gs[0, :]); lax.axis("off")
        lax.legend(handles=[
            Line2D([], [], marker="o", ls="", mfc="0.86", mec="0.55", ms=8, label="Absent"),
            Line2D([], [], marker="x", ls="", color="#1f1fff", mew=2, ms=8, label="Present"),
            Line2D([], [], marker="s", ls="", mfc="none", mec="darkorange", ms=8,
                   label="Masked SNP (not counted)")],
            title="T-to-C substitution at\na uridine site:", loc="center left",
            bbox_to_anchor=(0.15, 0.5), frameon=False, fontsize=9, title_fontsize=9)
        cax = lax.inset_axes([0.6, 0.35, 0.25, 0.18])
        cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP, norm=norm), cax=cax,
                          orientation="horizontal")
        cb.set_ticks(list(range(1, int(a.vmax) + 1)) if a.vmax <= 6 else
                     [t for t in (4, 8, 12, 16) if t <= a.vmax])
        cax.set_title("Coverage of T-to-C\nsubstitutions per uridine site:", fontsize=9)
        a.show_masked = True
        if not a.xlim:   # one common window for all panels (from both BAMs)
            allp = np.array([p for u in order for p in cov_all[u]]) + 1
            if a.full_range:
                lo, hi = allp.min(), allp.max()
            else:
                lo, hi = np.percentile(allp, [a.range_pct, 100 - a.range_pct])
            pad = max(15, int(0.03 * (hi - lo)))
            a.xlim = [int(lo) - pad, int(hi) + pad]
        it = f"$\\it{{{g}}}$"
        if a.compare_panels == "masked":
            panels = [(0, 0, cov_p, conv_p, snps, f"{it}: pre-consensus reads, with SNP mask"),
                      (0, 1, cov_c, conv_c, snps, f"{it}: consensus.bam, with SNP mask")]
        elif a.compare_panels == "nomask":
            panels = [(0, 0, cov_p, conv_p, set(), f"{it}: pre-consensus reads, no mask"),
                      (0, 1, cov_c, conv_c, set(), f"{it}: consensus.bam, no mask")]
        else:
            panels = [(0, 0, cov_p, conv_p, set(), f"{it}: pre-consensus reads, no mask"),
                      (0, 1, cov_p, conv_p, snps, f"{it}: pre-consensus reads, with SNP mask"),
                      (1, 0, cov_c, conv_c, set(), f"{it}: consensus.bam, no mask"),
                      (1, 1, cov_c, conv_c, snps, f"{it}: consensus.bam, with SNP mask")]
        print(f"{'panel':42s} transcripts  labeled  T>C events removed by mask")
        if a.reads_in_rows:   # every read of every UMI, per BAM
            reads_pre = collect_cell_reads(bam, fa, locs[g], cells[g], a)
            bq = a.min_bq
            a.min_bq = a.consensus_min_bq
            reads_con = collect_cell_reads(bam_c, fa, locs[g], cells[g], a)
            a.min_bq = bq
            mx = max((len(v) for v in reads_pre.values()), default=1)
            print(f"--reads-in-rows: up to {mx} pre-consensus reads per transcript")
        for r_, c2, cv, cn, m, title in panels:
            if a.reads_in_rows:
                a._cell_reads = {g: reads_pre if cv is cov_p else reads_con}
            ax = fig.add_subplot(gs[1 + r_, c2])
            o = [u for u in order]
            draw_panel(ax, g, "", locs[g], cv, cn, m, a, norm, order=o, title=title, mark=snps)
            if a.label_umis:
                ax.set_yticks(range(1, len(o) + 1))
                ax.set_yticklabels(o, fontsize=a.umi_fontsize, family="monospace")
                for t, u in zip(ax.get_yticklabels(), o):
                    if u in reps:
                        t.set_color("#e8302a"); t.set_fontweight("bold")
            from matplotlib.patches import Rectangle
            xl = ax.get_xlim()
            for u in reps:
                ri = o.index(u) + 1
                ax.axhspan(ri - 0.5, ri + 0.5, color="#e8302a", alpha=0.08, zorder=0)
                ax.add_patch(Rectangle((xl[0], ri - 0.5), xl[1] - xl[0], 1.0, fill=False,
                                       ls="--", ec="#e8302a", lw=1.3, zorder=4, clip_on=False))
                if not a.label_umis and c2 == 1:   # UMI name right of the right-hand panels
                    ax.annotate(u, xy=(xl[1], ri), xytext=(4, 0), textcoords="offset points",
                                ha="left", va="center", fontsize=7, color="#e8302a",
                                family="monospace", fontweight="bold", annotation_clip=False)
            shown = [u for u in o if cv.get(u)]
            nl = sum(1 for u in shown if any(k >= a.min_conv_reads and p not in m
                                              for p, k in cn[u].items()))
            rm = sum(1 for u in shown for p in cn[u] if p in m)
            name = title.split(": ", 1)[1]
            ax.set_ylabel(f"{len(shown)} transcripts ({nl} labeled, {len(shown) - nl} unlabeled)",
                          fontsize=9)
            print(f"{name:42s} {len(shown):11d}  {nl:7d}  {rm:d}")
        print("Rows are the same transcripts in all panels; an empty row = UMI not in that BAM.")
        label_cells(fig, [cells[g]])
        a.out = out_with_cell(a.out, [cells[g]])
        fig.savefig(a.out, dpi=300, bbox_inches="tight")
        print(f"Saved: {a.out}")
        return

    if a.paper:
        if len(a.genes) != 1:
            sys.exit("ERROR: --scnt-layout works with one gene")
        g = a.genes[0]
        cov, conv, snps = data[g]
        cr = collect_cell_reads(bam, fa, locs[g], cells[g], a)
        if a.xlim:
            lo, hi = a.xlim[0] - 1, a.xlim[1]
            cr = {u: [(r0, r1, {p for p in uri if lo <= p < hi}, {p for p in tc if lo <= p < hi})
                      for r0, r1, uri, tc in rs] for u, rs in cr.items()}
        region = "3'UTR" if a.utr3 else ("region" if a.xlim else "")
        if a.compare_layout == "mask" or (a.compare_consensus and a.compare_layout == "side"):
            # the same figure twice, side by side: input = pre-consensus BAM / consensus.bam;
            # same cell, same rows, same highlighted UMI, same window
            mask_mode = a.compare_layout == "mask"
            c_, s_, e_, st_ = locs[g]
            if mask_mode:        # same reads, counted without | with the SNP mask
                if not snps:
                    print(f"NOTE: no masked SNP site in this window of {g}; both halves identical")
                cov_c, conv_c, cr_c = cov, conv, cr
                a.show_masked = True
            else:
                bam_c = pysam.AlignmentFile(a.compare_consensus, "rb")
                bq = a.min_bq
                a.min_bq = a.consensus_min_bq
                cov_c, conv_c, _ = collect(bam_c, fa, c_, s_, e_, st_, cells[g], a)
                cr_c = collect_cell_reads(bam_c, fa, locs[g], cells[g], a)
                a.min_bq = bq
            if a.xlim:
                lo, hi = a.xlim[0] - 1, a.xlim[1]
                cov_c = defaultdict(set, {u: {p for p in ps if lo <= p < hi}
                                          for u, ps in cov_c.items()})
                conv_c = defaultdict(lambda: defaultdict(int),
                                     {u: defaultdict(int, {p: k for p, k in cs.items()
                                                           if lo <= p < hi})
                                      for u, cs in conv_c.items()})
                cr_c = {u: [(r0, r1, {p for p in uri if lo <= p < hi},
                             {p for p in tc if lo <= p < hi}) for r0, r1, uri, tc in rs]
                        for u, rs in cr_c.items()}
            order, lab, _ = order_umis(cov, conv, set() if mask_mode else snps,
                                       a.min_conv_reads, a.max_umis)
            if a.umi_reads and a.umi_reads != "auto":
                rep = a.umi_reads
            elif mask_mode and snps:   # highlight a transcript carrying a masked variant
                cand = [u for u in order if cr.get(u) and any(p in snps for p in conv[u])]
                cand = cand or [u for u in order if cr.get(u)]
                rep = max(cand, key=lambda u: (sum(p in snps for p in conv[u]), len(cr[u])))
            else:
                cand = [u for u in order if cr.get(u)]
                rep = max(cand, key=lambda u: (u in lab, len(cr[u]))) if cand else None
            if rep and rep not in order:
                order = [rep] + order[:-1]
            allp = np.array([p for u in order for p in (cov[u] | cov_c.get(u, set()))]) + 1
            if a.xlim:
                xl = a.xlim
            else:
                lo_, hi_ = np.percentile(allp, [a.range_pct, 100 - a.range_pct])
                pad = max(10, int(0.03 * (hi_ - lo_)))
                xl = [int(lo_) - pad, int(hi_) + pad]
            nrow = len(order)
            main_h = max(2.5, (a.row_height or 0.13) * nrow + 0.6)
            kmax = max(len(cr.get(rep, [])), len(cr_c.get(rep, [])), 1)
            inset_h = max(2.6, 0.22 * kmax + 0.9)
            a._paper_h = (main_h, inset_h)
            fig = plt.figure(figsize=(2 * a.width + 0.6, 0.9 + main_h + inset_h + 0.6))
            sf = fig.subfigures(1, 2, wspace=0.06)
            print(f"Highlighted transcript: UMI {rep}")
            base = a.sample_label + " — " if a.sample_label else ""
            if mask_mode:
                t1, t2, s1 = base + "counted without SNP mask", base + "counted with SNP mask", set()
            else:
                t1, t2, s1 = base + "pre-consensus BAM", base + "consensus.bam", snps
            plot_paper(a, g, locs[g], cells[g], cov, conv, s1, cr, region, fig=sf[0],
                       rep_override=rep, order_override=order, xlim_override=xl,
                       panel_title=t1)
            plot_paper(a, g, locs[g], cells[g], cov_c, conv_c, snps, cr_c, region, fig=sf[1],
                       rep_override=rep, order_override=order, xlim_override=xl,
                       panel_title=t2)
            label_cells(fig, [cells[g]])
            a.out = out_with_cell(a.out, [cells[g]])
            fig.savefig(a.out, dpi=300, bbox_inches="tight")
            print(f"Saved: {a.out}")
            return
        crc = None
        if a.compare_consensus:   # consensus record(s) of the representative UMI in the inset
            bq = a.min_bq
            a.min_bq = a.consensus_min_bq
            crc = collect_cell_reads(pysam.AlignmentFile(a.compare_consensus, "rb"), fa,
                                     locs[g], cells[g], a)
            a.min_bq = bq
            if a.xlim:
                crc = {u: [(r0, r1, {p for p in uri if lo <= p < hi},
                            {p for p in tc if lo <= p < hi}) for r0, r1, uri, tc in rs]
                       for u, rs in crc.items()}
        plot_paper(a, g, locs[g], cells[g], cov, conv, snps, cr, region, cons_reads=crc)
        return

    norm = Normalize(vmin=1, vmax=a.vmax)
    if a.before_after:
        if len(a.genes) != 1 or not a.snp_csv:
            sys.exit("ERROR: --before-after needs exactly one gene and --snp-csv")
        g = a.genes[0]
        cov, conv, snps = data[g]
        a.show_masked = True
        # (mask used for counting, title, mask positions drawn as lines, short name)
        panels = [(set(), f"$\\it{{{g}}}$: counted without SNP mask", snps, "no mask")]
        c_, s_, e_, _ = locs[g]
        if a.old_snp_csv and not a.control_mask:
            a.control_mask, a.control_label = a.old_snp_csv, "OLD 1-based mask, as dynast reads it"
        for k in a.control_shift:
            sh = {p + k for p in snps}
            panels.append((sh, f"$\\it{{{g}}}$: negative control - mask shifted by {k:+d} bp",
                           sh, f"mask shifted {k:+d} bp"))
        if a.control_mask:
            cm = {p for p in load_snps(a.control_mask, c_, "0based") if s_ <= p < e_}
            panels.append((cm, f"$\\it{{{g}}}$: negative control - {a.control_label}", cm,
                           a.control_label))
        panels.append((snps, f"$\\it{{{g}}}$: counted with WGS SNP mask (0-based)", snps,
                       "new 0-based mask"))
        if a.only_masked:
            panels = panels[-1:]
        order, _, _ = order_umis(cov, conv, set(), a.min_conv_reads, a.max_umis)
        rowh = a.row_height or (0.22 if a.reads_in_rows else 0.11)
        h = max(2.0, rowh * len(order) + 1.0)
        heights = [h] * len(panels)
    else:
        rowh = a.row_height or (0.22 if a.reads_in_rows else 0.11)
        heights = [max(2.0, rowh * min(len(data[g][0]), a.max_umis) + 1.0) for g in a.genes]
    umi_reads = None
    if a.umi_reads:
        g0 = a.genes[0]
        cov0, conv0, snps0 = data[g0]
        cell0 = cells[g0] if not a.before_after else norm_bc(a.cell)
        if a.umi_reads == "auto":
            c_, s_, e_, st_ = locs[g0]
            n_by = defaultdict(int)
            for r in iter_gene_reads(bam, c_, s_, e_, st_, a.strand, a.cb_tag, a.ub_tag,
                                     a.min_mapq):
                if norm_bc(r.get_tag(a.cb_tag)) == cell0:
                    n_by[r.get_tag(a.ub_tag)] += 1
            labeled0 = {u for u in conv0 if any(p not in snps0 for p in conv0[u])}
            pick = sorted(n_by, key=lambda u: (u in labeled0, n_by[u]), reverse=True)
            if not pick:
                sys.exit("ERROR: no reads for --umi-reads")
            umi0 = pick[0]
            print(f"UMIs with most reads in this cell: "
                  + ", ".join(f"{u}={n_by[u]}" for u in sorted(n_by, key=n_by.get, reverse=True)[:5]))
        else:
            umi0 = a.umi_reads
        umi_reads = collect_umi_reads(bam, fa, locs[g0], cell0, umi0, a, snps0)
        if not umi_reads:
            sys.exit(f"ERROR: UMI {umi0} has no reads for {g0} in cell {cell0}")
        if all(not any(b for b in r[0]) for r in umi_reads):
            sys.exit("ERROR: reads without alignment blocks")
        print(f"--umi-reads: UMI {umi0} in cell {cell0}: {len(umi_reads)} reads"
              + (" (consensus BAM: usually 1 record per UMI -- use the pre-consensus BAM)"
                 if len(umi_reads) <= 2 else ""))
        rh = max(1.2, 0.18 * len(umi_reads) + 0.8)
        heights = [rh] if a.umi_only else heights + [rh]
    fig_h = sum(heights) + 0.6 + 1.0 * len(heights) + 0.9
    fig = plt.figure(figsize=(a.width, fig_h))
    gs = fig.add_gridspec(len(heights) + 1, 1, height_ratios=[0.6] + heights,
                          hspace=1.0 / np.mean([0.6] + heights),   # ~1 inch between panels
                          top=1 - 0.3 / fig_h, bottom=0.6 / fig_h)  # fixed margins (inches)

    # legend row (as in the paper): markers on the left, colour bar on the right
    lax = fig.add_subplot(gs[0]); lax.axis("off")
    handles = [
        Line2D([], [], marker="o", ls="", mfc="0.86", mec="0.55", ms=8, label="Absent"),
        Line2D([], [], marker="x", ls="", color="#1f1fff", mew=2, ms=8, label="Present")]
    if a.show_masked or a.before_after:
        handles.append(Line2D([], [], marker="s", ls="", mfc="none", mec="darkorange", ms=8,
                              label="Masked SNP (not counted)"))
    lax.legend(handles=handles, title="T-to-C substitution at\na uridine site:", loc="center left",
               bbox_to_anchor=(0.12, 0.5), frameon=False, fontsize=9, title_fontsize=9)
    cax = lax.inset_axes([0.6, 0.35, 0.3, 0.18])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP, norm=norm), cax=cax, orientation="horizontal")
    cb.set_ticks(list(range(1, int(a.vmax) + 1)) if a.vmax <= 6 else
                 [t for t in (4, 8, 12, 16, 20, 24) if t <= a.vmax])
    cax.set_title("Coverage of T-to-C\nsubstitutions per uridine site:", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    if umi_reads and a.umi_only:
        pass
    elif a.before_after:
        for i, (m, title, mark, what) in enumerate(panels):
            ax = fig.add_subplot(gs[i + 1])
            nl, nu = draw_panel(ax, g, "", locs[g], cov, conv, m, a, norm,
                                order=order, title=title, mark=mark)
            tot = nl + nu
            n_rm = sum(1 for u in order for p in conv[u] if p in m)
            print(f"{g} {what:36s}: {tot} UMIs, {nl} labeled ({(nl / tot if tot else 0):.0%}), "
                  f"{n_rm} T>C events removed")
    else:
        for i, g in enumerate(a.genes):
            ax = fig.add_subplot(gs[i + 1])
            cov, conv, snps = data[g]
            nl, nu = draw_panel(ax, g, labels[i], locs[g], cov, conv, snps, a, norm)
            if len(set(cells.values())) > 1:
                ax.set_title(ax.get_title() + f"\ncell {cells[g]}", fontsize=10)
            tot = nl + nu
            print(f"{g}: {tot} UMIs, {nl} labeled ({(nl / tot if tot else 0):.0%}), {nu} unlabeled"
                  + (f", {len(snps)} masked SNP sites in gene" if a.snp_csv else ""))

    if umi_reads:
        ax = fig.add_subplot(gs[len(heights)])
        allp = [p for b, _, _ in umi_reads for blk in b for p in blk]
        pad = max(10, int(0.05 * (max(allp) - min(allp))))
        draw_umi_reads(ax, umi_reads, locs[a.genes[0]][0], umi0, a,
                       xlim=(min(allp) + 1 - pad, max(allp) + pad))

    label_cells(fig, [cells[g] for g in a.genes])
    a.out = out_with_cell(a.out, [cells[g] for g in a.genes])
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    print(f"Saved: {a.out}")


if __name__ == "__main__":
    main()
