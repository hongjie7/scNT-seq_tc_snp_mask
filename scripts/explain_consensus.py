#!/usr/bin/env python3
"""
explain_consensus.py -- show how dynast consensus turns the reads of ONE molecule
(cell barcode + UMI) into one consensus read, using the two real BAM files.

Panel A  the reads of the UMI in the pre-consensus BAM and its record in
         consensus.bam, along the genome; every mismatch to the reference is
         marked (T>C on the counting strand as X, other mismatches as triangles)
Panel B  a table of every position where any read (or the consensus) differs
         from the reference: base + base quality of each read, the quality-
         weighted vote, the base dynast's rule predicts, the base actually in
         consensus.bam, and whether dynast count would count it as T>C

dynast consensus rule (dynast 1.0.1, preprocessing/consensus.py):
  at each position, add up the base qualities of all reads per base; the base
  with the highest sum wins; if that sum is < 27 (or tied) the reference base
  is used; consensus quality = min(sum, 42).
dynast count then counts a T>C only if the consensus quality is > 27.

Usage:
  python explain_consensus.py --bam split_bams/<donor>.bam \
      --consensus consensus.bam --fasta genome.fa --gtf annotation.gtf \
      --gene SOX2 [--cell CB] [--umi UB] --out explain_SOX2.pdf
Without --cell/--umi, the UMI with the most reads that shows a disagreement
between its reads (or with the consensus) is chosen.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_tc_panel import find_genes  # noqa: E402

QMIN, QCAP = 27, 42
BLUE, ORANGE, GREY = "#1f1fff", "#d97706", "0.55"


def norm_bc(bc):
    return bc[:-2] if bc.endswith("-1") else bc


def read_bases(r):
    """{ref_pos0: (base, qual)} for aligned (match/mismatch) positions."""
    seq, qual = r.query_sequence, r.query_qualities
    return {rp: (seq[qp], qual[qp] if qual is not None else 30)
            for qp, rp in r.get_aligned_pairs(matches_only=True)}


def group_reads(bam, chrom, start, end, cb_tag, ub_tag, mapq, cell=None):
    """{(cell, umi): [read, ...]} for reads in the region."""
    g = defaultdict(list)
    for r in bam.fetch(chrom, start, end):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < mapq:
            continue
        if not (r.has_tag(cb_tag) and r.has_tag(ub_tag)):
            continue
        cb, ub = norm_bc(r.get_tag(cb_tag)), r.get_tag(ub_tag)
        if cb in ("-", "") or ub in ("-", ""):
            continue
        if cell and cb != cell:
            continue
        g[(cb, ub)].append(r)
    return g


def predict(bases_list, refb):
    """dynast rule for one position: (consensus base, summed quality, votes)."""
    votes = defaultdict(int)
    for b, q in bases_list:
        if b != "N":
            votes[b] += q
    if not votes:
        return refb, 0, votes
    best = max(votes.values())
    winners = [b for b, v in votes.items() if v == best]
    if best < QMIN or len(winners) > 1:
        return refb, min(votes.get(refb, 0), QCAP), votes
    return winners[0], min(best, QCAP), votes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="pre-consensus BAM (indexed)")
    ap.add_argument("--consensus", required=True, help="dynast consensus.bam (indexed)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", help="GTF, needed with --gene")
    ap.add_argument("--gene", help="gene name (region and strand from the GTF)")
    ap.add_argument("--region", help="chr:start-end instead of --gene (1-based)")
    ap.add_argument("--strand", choices=["+", "-"], help="gene strand with --region")
    ap.add_argument("--cell", help="cell barcode")
    ap.add_argument("--umi", help="UMI")
    ap.add_argument("--min-reads", type=int, default=3,
                    help="auto choice: UMI needs >= N reads (default 3)")
    ap.add_argument("--max-cols", type=int, default=30, help="max positions in the table")
    ap.add_argument("--cb-tag", default="CB")
    ap.add_argument("--ub-tag", default="UB")
    ap.add_argument("--min-mapq", type=int, default=255)
    ap.add_argument("--out", default="explain_consensus.pdf")
    a = ap.parse_args()

    pre = pysam.AlignmentFile(a.bam, "rb")
    con = pysam.AlignmentFile(a.consensus, "rb")
    fa = pysam.FastaFile(a.fasta)
    if a.gene:
        if not a.gtf:
            sys.exit("ERROR: --gene needs --gtf")
        chrom, start, end, strand = find_genes(a.gtf, [a.gene], set(pre.references))[a.gene]
        name = a.gene
    elif a.region:
        chrom, rng = a.region.replace(",", "").split(":")
        s0, e0 = rng.split("-")
        start, end, strand, name = int(s0) - 1, int(e0), a.strand or "+", a.region
    else:
        sys.exit("ERROR: give --gene or --region")
    ref_u, conv = ("T", "C") if strand == "+" else ("A", "G")

    groups = group_reads(pre, chrom, start, end, a.cb_tag, a.ub_tag, a.min_mapq, a.cell)
    cgroups = group_reads(con, chrom, start, end, a.cb_tag, a.ub_tag, a.min_mapq, a.cell)

    def cons_record(key):
        rs = cgroups.get(key, [])
        tagged = [r for r in rs if r.has_tag("RN")]
        return (tagged or rs)[0] if rs else None

    def disagreement(key):
        rs = groups[key]
        per = defaultdict(list)
        for r in rs:
            for p, bq in read_bases(r).items():
                per[p].append(bq[0])
        n = sum(1 for p, bs in per.items() if len(set(bs)) > 1)
        cr = cons_record(key)
        if cr is not None:
            refseq = fa.fetch(chrom, cr.reference_start, cr.reference_end).upper()
            n += sum(1 for p, (b, q) in read_bases(cr).items()
                     if b != refseq[p - cr.reference_start])
        return n

    if a.umi:
        if not a.cell:
            keys = [k for k in groups if k[1] == a.umi]
            if not keys:
                sys.exit(f"ERROR: UMI {a.umi} not found in {name}")
            key = max(keys, key=lambda k: len(groups[k]))
        else:
            key = (norm_bc(a.cell), a.umi)
            if key not in groups:
                sys.exit(f"ERROR: cell {a.cell} UMI {a.umi} has no reads in {name}")
    else:
        cand = [k for k in groups if len(groups[k]) >= a.min_reads and cons_record(k) is not None]
        if not cand:
            sys.exit(f"ERROR: no UMI with >= {a.min_reads} reads and a consensus record in "
                     f"{name}; lower --min-reads")
        key = max(cand, key=lambda k: (disagreement(k), len(groups[k])))
    cell, umi = key
    reads = sorted(groups[key], key=lambda r: (r.reference_start, r.reference_end))
    crec = cons_record(key)
    print(f"{name} ({chrom}, {strand} strand; uridine = genomic {ref_u}, T>C = genomic "
          f"{ref_u}>{conv})")
    print(f"cell {cell}  UMI {umi}: {len(reads)} pre-consensus reads, consensus record: "
          + ("none" if crec is None else
             f"RN={crec.get_tag('RN') if crec.has_tag('RN') else 'none (single read)'}"))

    rb = [read_bases(r) for r in reads]
    cb_ = read_bases(crec) if crec is not None else {}
    lo = min(r.reference_start for r in reads + ([crec] if crec else []))
    hi = max(r.reference_end for r in reads + ([crec] if crec else []))
    refseq = fa.fetch(chrom, lo, hi).upper()
    ref = lambda p: refseq[p - lo]

    # informative positions: any read or the consensus differs from the reference
    inf = sorted({p for d in rb for p, (b, q) in d.items() if b != ref(p) and b != "N"} |
                 {p for p, (b, q) in cb_.items() if b != ref(p)})
    if len(inf) > a.max_cols:
        score = {p: sum(1 for d in rb if p in d and d[p][0] != ref(p)) for p in inf}
        inf = sorted(sorted(inf, key=lambda p: -score[p])[: a.max_cols])
        print(f"NOTE: showing the {a.max_cols} positions with the most mismatching reads")

    rows = []
    for p in inf:
        bl = [d[p] for d in rb if p in d]
        pb, pq, votes = predict(bl, ref(p))
        cbase, cq = cb_.get(p, ("·", None))
        is_tc = ref(p) == ref_u and cbase == conv
        counted = is_tc and cq is not None and cq > QMIN
        rows.append((p, bl, pb, pq, votes, cbase, cq, counted))

    # ---------------------------------------------------------------- figure
    ncol = max(len(inf), 1)
    nread = len(reads)
    hA = 0.32 * (nread + 1) + 1.0
    hB = 0.3 * (nread + 4) + 1.2
    width = max(8.0, 0.55 * ncol + 2.6)
    fig = plt.figure(figsize=(width, hA + hB + 1.4))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.7, hA, hB], hspace=0.9 / np.mean([0.7, hA, hB]))

    hax = fig.add_subplot(gs[0]); hax.axis("off")
    hax.text(0.5, 1.0, f"How dynast consensus merges the reads of one molecule — "
             f"{name}, cell {cell}, UMI {umi}", ha="center", va="top", fontsize=11,
             fontweight="bold", transform=hax.transAxes)
    hax.text(0.5, 0.45, "Per position, base qualities are summed per base; the highest sum wins; "
             f"if it is < {QMIN} (or tied) the reference base is used; consensus quality = "
             f"min(sum, {QCAP}). dynast count counts a T>C only if quality > {QMIN}.",
             ha="center", va="top", fontsize=8.5, color="0.25", wrap=True,
             transform=hax.transAxes)

    # Panel A
    ax = fig.add_subplot(gs[1])
    ys = list(range(1, nread + 1)) + [nread + 1.6]
    for j, (r, d) in enumerate(zip(reads + ([crec] if crec else []), rb + ([cb_] if crec else []))):
        y = ys[j]
        is_c = crec is not None and j == nread
        for b0, b1 in r.get_blocks():
            ax.plot([b0 + 1, b1], [y, y], color=ORANGE if is_c else "0.5", lw=2.2 if is_c else 1.4,
                    solid_capstyle="butt", zorder=1)
        for p, (b, q) in d.items():
            if b == ref(p) or b == "N":
                continue
            tc = ref(p) == ref_u and b == conv
            low = q is not None and q <= QMIN
            ax.scatter([p + 1], [y], marker="x" if tc else "^", s=46 if tc else 26,
                       color=BLUE if tc else "#c2410c", linewidths=1.6, zorder=3,
                       alpha=0.4 if low else 1.0)
    for k, p in enumerate(inf, 1):
        ax.axvline(p + 1, color="0.75", lw=0.6, ls=":", zorder=0)
        ax.text(p + 1, 0.15, str(k), ha="center", va="bottom", fontsize=6.5, color="0.35")
    ax.set_ylim(ys[-1] + 0.8, -0.2)
    ax.set_yticks(ys[: nread + (1 if crec else 0)])
    ax.set_yticklabels([f"read {i}" for i in range(1, nread + 1)] + (["consensus"] if crec else []),
                       fontsize=8)
    if crec:
        ax.get_yticklabels()[-1].set_color(ORANGE)
        ax.get_yticklabels()[-1].set_fontweight("bold")
        ax.axhline(nread + 0.8, color="0.4", lw=0.8)
    ax.set_xlim(lo - 3, hi + 3)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.set_xlabel(f"{chrom} position", fontsize=9)
    ax.set_title("A   Reads in the pre-consensus BAM (grey) and the record in consensus.bam "
                 "(orange); numbers = columns in table B", fontsize=9.5, loc="left")
    ax.legend(handles=[
        Line2D([], [], marker="x", ls="", color=BLUE, mew=1.6, label=f"T>C ({ref_u}>{conv} on genome)"),
        Line2D([], [], marker="^", ls="", color="#c2410c", label="other mismatch"),
        Line2D([], [], marker="x", ls="", color=BLUE, alpha=0.4, mew=1.6,
               label=f"faded: base quality ≤ {QMIN}")],
        loc="upper right", fontsize=7.5, frameon=True, framealpha=0.9)

    # Panel B: table
    bx = fig.add_subplot(gs[2]); bx.axis("off")
    bx.set_title("B   Every position where a read or the consensus differs from the reference",
                 fontsize=9.5, loc="left")
    labels = ([f"read {i}" for i in range(1, nread + 1)] +
              ["vote (sum of Q)", "dynast rule →", "consensus.bam", "counted as T>C"])
    nr = len(labels)
    cw, rh, x0 = 1.0, 1.0, 2.6
    bx.set_xlim(0, x0 + ncol * cw + 0.2)
    bx.set_ylim(nr + 1.6, -0.2)
    for i, lab in enumerate(labels):
        bx.text(x0 - 0.15, i + 1.5, lab, ha="right", va="center", fontsize=7.5,
                color=ORANGE if lab == "consensus.bam" else "black",
                fontweight="bold" if i >= nread else "normal")
    for k, (p, bl, pb, pq, votes, cbase, cq, counted) in enumerate(rows):
        xc = x0 + k * cw
        bx.text(xc + 0.5, 0.25, f"{k + 1}\n{ref(p)}", ha="center", va="bottom", fontsize=6.5,
                color="0.3")
        for i, d in enumerate(rb):
            if p in d:
                b, q = d[p]
                mis = b != ref(p)
                tc = ref(p) == ref_u and b == conv
                fc = ("#dbe4ff" if tc else "#ffe4cc") if mis else "0.96"
                bx.add_patch(Rectangle((xc, i + 1), cw, rh, fc=fc, ec="white", lw=0.8))
                bx.text(xc + 0.5, i + 1.5, f"{b}{q}", ha="center", va="center", fontsize=6.3,
                        color=(BLUE if tc else "#c2410c") if mis else "0.45",
                        alpha=0.45 if q <= QMIN else 1.0)
            else:
                bx.add_patch(Rectangle((xc, i + 1), cw, rh, fc="white", ec="0.92", lw=0.5))
        y = nread + 1
        vt = "\n".join(f"{b}{v}" for b, v in sorted(votes.items(), key=lambda kv: -kv[1])[:2])
        bx.add_patch(Rectangle((xc, y), cw, rh, fc="0.93", ec="white"))
        bx.text(xc + 0.5, y + 0.5, vt, ha="center", va="center", fontsize=5.6)
        y += 1
        pmis = pb != ref(p)
        bx.add_patch(Rectangle((xc, y), cw, rh, fc="#e5e7eb", ec="white"))
        bx.text(xc + 0.5, y + 0.5, f"{pb}", ha="center", va="center", fontsize=7,
                color=BLUE if pmis else "0.3", fontweight="bold" if pmis else "normal")
        y += 1
        agree = cbase == pb or cbase == "·"
        bx.add_patch(Rectangle((xc, y), cw, rh, fc="#fff1e0" if agree else "#fecaca", ec="white"))
        bx.text(xc + 0.5, y + 0.5, f"{cbase}{cq if cq is not None else ''}", ha="center",
                va="center", fontsize=6.3, color=ORANGE, fontweight="bold")
        y += 1
        bx.add_patch(Rectangle((xc, y), cw, rh, fc="#dbe4ff" if counted else "white", ec="0.9"))
        bx.text(xc + 0.5, y + 0.5, "yes" if counted else "–", ha="center", va="center",
                fontsize=6.5, color=BLUE if counted else "0.5", fontweight="bold" if counted else "normal")
    bx.text(0.0, nr + 1.35, "Cells: base + base quality. Blue = T>C on the counting strand, orange = "
            "other mismatch, grey = reference; faded = quality ≤ 27.\n'·' = position not in the "
            "consensus record. Red consensus cell = differs from the rule's prediction.",
            fontsize=6.8, color="0.3", va="top")

    # terminal summary
    print(f"\n{'col':>3} {'position':>14} ref  reads (base+Q)                      vote          "
          "rule  consensus  T>C counted")
    for k, (p, bl, pb, pq, votes, cbase, cq, counted) in enumerate(rows, 1):
        rs = " ".join(f"{b}{q}" for b, q in bl)
        vt = " ".join(f"{b}:{v}" for b, v in sorted(votes.items(), key=lambda kv: -kv[1]))
        print(f"{k:>3} {chrom}:{p + 1:>10,}  {ref(p)}   {rs[:34]:34s} {vt[:13]:13s}  {pb:>3}  "
              f"{cbase + (str(cq) if cq is not None else ''):>9}  {'yes' if counted else '-'}")
    mism = sum(1 for (p, bl, pb, pq, v, cbase, cq, c) in rows if cbase not in (pb, "·"))
    if mism:
        print(f"NOTE: {mism} position(s) where consensus.bam differs from the rule's prediction "
              "(e.g. different quality filters or alignment); shown in red")

    root, ext = os.path.splitext(a.out)
    out = f"{root}_{cell}_{umi}{ext}"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
