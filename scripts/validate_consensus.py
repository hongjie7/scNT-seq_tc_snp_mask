#!/usr/bin/env python3
"""
validate_consensus.py -- check with the two real BAM files, over MANY molecules,
what dynast consensus did with T>C seen in the pre-consensus reads.

For every molecule (cell barcode + UMI) of the given genes that has a record in
consensus.bam, and every uridine position covered by that record and by >= 1 of
its pre-consensus reads:
  n = pre-consensus reads covering the position
  k = of these, reads with T>C (genomic T>C on + strand genes, A>G on - strand)
  kept    = the consensus record has the converted base there
  counted = kept AND consensus base quality > 27 (what dynast count uses)
It also checks, at every covered position, whether consensus.bam agrees with
dynast's rule applied to the pre-consensus reads (quality-weighted vote, sum
>= 27, else reference).

Figure:
  A  fraction of T>C kept in the consensus, by n (reads covering) and k (reads
     with T>C); numbers = positions
  B  how often a T>C seen in the reads ends up counted, by n
  C  agreement of consensus.bam with dynast's rule, by n

Usage:
  python validate_consensus.py --bam split_bams/<donor>.bam \
      --consensus consensus.bam --fasta genome.fa --gtf annotation.gtf \
      --genes SOX2 ACTB GAPDH --out validate_consensus.pdf
"""
import argparse
import os
import sys
from collections import defaultdict, Counter

import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_tc_panel import find_genes  # noqa: E402

QMIN, QCAP = 27, 42


def norm_bc(bc):
    return bc[:-2] if bc.endswith("-1") else bc


def groups(bam, chrom, start, end, a):
    g = defaultdict(list)
    for r in bam.fetch(chrom, start, end):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < a.min_mapq:
            continue
        if not (r.has_tag(a.cb_tag) and r.has_tag(a.ub_tag)):
            continue
        cb, ub = norm_bc(r.get_tag(a.cb_tag)), r.get_tag(a.ub_tag)
        if cb in ("-", "") or ub in ("-", ""):
            continue
        g[(cb, ub)].append(r)
    return g


def bases(r):
    s, q = r.query_sequence, r.query_qualities
    return {rp: (s[qp], q[qp] if q is not None else 30)
            for qp, rp in r.get_aligned_pairs(matches_only=True)}


def rule(bl, refb):
    v = defaultdict(int)
    for b, q in bl:
        if b != "N":
            v[b] += q
    if not v:
        return refb
    best = max(v.values())
    w = [b for b, x in v.items() if x == best]
    return refb if best < QMIN or len(w) > 1 else w[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="pre-consensus BAM")
    ap.add_argument("--consensus", required=True, help="dynast consensus.bam")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--max-n", type=int, default=6, help="pool n >= this into one column")
    ap.add_argument("--max-umis", type=int, default=0, help="per gene, 0 = all")
    ap.add_argument("--cb-tag", default="CB")
    ap.add_argument("--ub-tag", default="UB")
    ap.add_argument("--min-mapq", type=int, default=255)
    ap.add_argument("--out", default="validate_consensus.pdf")
    ap.add_argument("--table", help="optional TSV with one line per T>C position")
    a = ap.parse_args()

    pre = pysam.AlignmentFile(a.bam, "rb")
    con = pysam.AlignmentFile(a.consensus, "rb")
    fa = pysam.FastaFile(a.fasta)
    locs = find_genes(a.gtf, a.genes, set(pre.references))

    tc = defaultdict(lambda: [0, 0, 0])          # (n, k) -> [positions, kept, counted]
    agree = defaultdict(lambda: [0, 0])          # n -> [positions, agree with rule]
    nmol, molreads = Counter(), defaultdict(list)
    rows = []
    for g in a.genes:
        chrom, start, end, strand = locs[g]
        ref_u, conv = ("T", "C") if strand == "+" else ("A", "G")
        gp, gc = groups(pre, chrom, start, end, a), groups(con, chrom, start, end, a)
        keys = [k for k in gp if k in gc]
        if a.max_umis:
            keys = sorted(keys, key=lambda k: -len(gp[k]))[: a.max_umis]
        print(f"{g}: {len(gp):,} molecules in pre-consensus BAM, {len(gc):,} in consensus.bam, "
              f"{len(keys):,} in both")
        for key in keys:
            reads = gp[key]
            crs = gc[key]
            crec = ([r for r in crs if r.has_tag("RN")] or crs)[0]
            cb = bases(crec)
            rb = [bases(r) for r in reads]
            lo, hi = crec.reference_start, crec.reference_end
            refseq = fa.fetch(chrom, lo, hi).upper()
            nmol[g] += 1
            molreads[g].append(len(reads))
            for p, (cbase, cq) in cb.items():
                if not (start <= p < end):
                    continue
                refb = refseq[p - lo]
                bl = [d[p] for d in rb if p in d]
                n = len(bl)
                if n == 0:
                    continue
                nn = min(n, a.max_n)
                agree[nn][0] += 1
                agree[nn][1] += cbase == rule(bl, refb)
                if refb != ref_u:
                    continue
                k = sum(1 for b, q in bl if b == conv)
                if k == 0:
                    continue
                kept = cbase == conv
                counted = kept and cq > QMIN
                kk = min(k, nn)
                tc[(nn, kk)][0] += 1
                tc[(nn, kk)][1] += kept
                tc[(nn, kk)][2] += counted
                rows.append((g, key[0], key[1], chrom, p + 1, n, k, int(kept), int(counted)))

    if not rows:
        sys.exit("No T>C seen in pre-consensus reads at positions covered by a consensus record.")

    # ---------------------------------------------------------------- summary
    tot = sum(v[0] for v in tc.values())
    print(f"\nT>C positions in pre-consensus reads (covered by the consensus record): {tot:,}")
    print(f"{'n reads':>8} {'k with T>C':>11} {'positions':>10} {'kept':>7} {'counted':>8}")
    for (n, k) in sorted(tc):
        v = tc[(n, k)]
        nl = f"{n}+" if n == a.max_n else str(n)
        kl = f"{k}+" if k == a.max_n else str(k)
        print(f"{nl:>8} {kl:>11} {v[0]:>10,} {v[1] / v[0]:>7.0%} {v[2] / v[0]:>8.0%}")
    ag = sum(v[1] for v in agree.values()) / max(sum(v[0] for v in agree.values()), 1)
    print(f"\nconsensus.bam agrees with dynast's rule at {ag:.2%} of covered positions")

    # ---------------------------------------------------------------- figure
    N = a.max_n
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8),
                            gridspec_kw={"width_ratios": [1.25, 1, 1], "wspace": 0.35})
    ax = axs[0]
    mat = np.full((N, N), np.nan)
    for (n, k), v in tc.items():
        mat[k - 1, n - 1] = v[1] / v[0]
    im = ax.imshow(mat, origin="lower", cmap="RdBu_r", vmin=0, vmax=1, aspect="auto")
    for (n, k), v in tc.items():
        ax.text(n - 1, k - 1, f"{v[1] / v[0]:.0%}\n(n={v[0]:,})", ha="center", va="center",
                fontsize=7, color="white" if abs(v[1] / v[0] - 0.5) > 0.3 else "black")
    lab = [str(i) for i in range(1, N)] + [f"{N}+"]
    ax.set_xticks(range(N)); ax.set_xticklabels(lab)
    ax.set_yticks(range(N)); ax.set_yticklabels(lab)
    ax.set_xlabel("pre-consensus reads covering the position")
    ax.set_ylabel("of these, reads with T>C")
    ax.set_title("A  T>C kept in consensus.bam", fontsize=11, loc="left")
    fig.colorbar(im, ax=ax, fraction=0.046, label="fraction kept")

    ax = axs[1]
    ns = sorted({n for n, _ in tc})
    single = [sum(v[2] for (n, k), v in tc.items() if n == m and k == 1) /
              max(sum(v[0] for (n, k), v in tc.items() if n == m and k == 1), 1) for m in ns]
    allk = [sum(v[2] for (n, k), v in tc.items() if n == m and k == n) /
            max(sum(v[0] for (n, k), v in tc.items() if n == m and k == n), 1) for m in ns]
    x = np.arange(len(ns))
    ax.bar(x - 0.2, single, 0.4, color="#9ca3af", label="T>C in 1 read only")
    ax.bar(x + 0.2, allk, 0.4, color="#1f1fff", label="T>C in every covering read")
    ax.set_xticks(x); ax.set_xticklabels([f"{m}+" if m == N else str(m) for m in ns])
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("pre-consensus reads covering the position")
    ax.set_ylabel("fraction counted by dynast (kept, Q > 27)")
    ax.set_title("B  Single-read T>C vs T>C in all reads", fontsize=11, loc="left")
    ax.legend(fontsize=8, loc="lower left")

    ax = axs[2]
    na = sorted(agree)
    ax.bar(range(len(na)), [agree[m][1] / agree[m][0] for m in na], color="#d97706")
    for i, m in enumerate(na):
        ax.text(i, agree[m][1] / agree[m][0] + 0.005, f"{agree[m][0]:,}", ha="center",
                fontsize=7, rotation=90, va="bottom")
    ax.set_xticks(range(len(na))); ax.set_xticklabels([f"{m}+" if m == N else str(m) for m in na])
    lo_ = min(agree[m][1] / agree[m][0] for m in na)
    ax.set_ylim(max(0, lo_ - 0.05), 1.03)
    ax.set_xlabel("pre-consensus reads covering the position")
    ax.set_ylabel("consensus base = dynast rule")
    ax.set_title(f"C  consensus.bam vs rule ({ag:.2%} agree)", fontsize=11, loc="left")

    mol = ", ".join(f"{g}: {nmol[g]:,} molecules, mean {np.mean(molreads[g]):.1f} reads"
                    for g in a.genes if nmol[g])
    fig.suptitle(f"dynast consensus checked on the two BAM files ({mol})", fontsize=11, y=1.03)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    print(f"Saved: {a.out}")

    if a.table:
        with open(a.table, "w") as f:
            f.write("gene\tcell\tumi\tchrom\tpos_1based\tn_reads\tk_TC\tkept\tcounted\n")
            for r in rows:
                f.write("\t".join(map(str, r)) + "\n")
        print(f"Saved: {a.table}")


if __name__ == "__main__":
    main()
