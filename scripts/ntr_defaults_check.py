#!/usr/bin/env python3
"""
ntr_defaults_check.py -- look at your dynast count adata.h5ad and suggest
reasonable values for the NTR options of plot_tc_panel.py.

Prints:
  1. dataset size and overall NTR (labeled / all UMIs)
  2. how many genes pass different expression filters, and their NTR spread
  3. labeled UMIs of high-NTR genes (for --min-gene-labeled)
  4. UMIs per cell of well-expressed genes (for --min-cell-umis)
  5. suggested option values and a ready-to-run command

Usage:
  python ntr_defaults_check.py --adata 02_dynast_count/TC_rm_bg_0based/adata.h5ad
"""
import argparse
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--adata", required=True)
    ap.add_argument("--no-mt-ribo", action="store_true", help="skip MT-, RPL*, RPS* genes")
    a = ap.parse_args()

    import anndata
    ad = anndata.read_h5ad(a.adata)
    lab = unl = None
    for l1, l2 in (("labeled_TC", "unlabeled_TC"), ("X_l_TC", "X_n_TC")):
        if l1 in ad.layers and l2 in ad.layers:
            lab, unl = sp.csc_matrix(ad.layers[l1]), sp.csc_matrix(ad.layers[l2])
            break
    if lab is None:
        sys.exit("ERROR: no labeled/unlabeled TC layers")
    tot = lab + unl
    names = (ad.var["gene_name"].astype(str).values if "gene_name" in ad.var
             else ad.var_names.astype(str).values)
    df = pd.DataFrame({"gene": names,
                       "UMIs": np.asarray(tot.sum(0)).ravel(),
                       "labeled": np.asarray(lab.sum(0)).ravel(),
                       "cells": np.asarray((tot > 0).sum(0)).ravel()})
    df["NTR"] = np.where(df.UMIs > 0, df.labeled / df.UMIs.clip(lower=1), 0)
    if a.no_mt_ribo:
        df = df[~df.gene.str.match(r"^(MT-|RPL|RPS)")]
    ncell = ad.n_obs
    overall = df.labeled.sum() / df.UMIs.sum()
    print(f"1. {ncell:,} cells, {ad.n_vars:,} genes, {int(df.UMIs.sum()):,} UMIs; "
          f"overall NTR {overall:.2%}")

    min_cells = max(20, int(0.01 * ncell))
    print(f"\n2. Genes passing filters (expressed in >= {min_cells} cells = max(20, 1% of cells)):")
    print(f"   {'min UMIs':>9} {'genes':>7} {'NTR 10%':>8} {'median':>8} {'NTR 90%':>8}")
    rows = {}
    for m in (50, 100, 200, 500, 1000, 2000):
        p = df[(df.UMIs >= m) & (df.cells >= min_cells)]
        if len(p) == 0:
            continue
        q = p.NTR.quantile([0.1, 0.5, 0.9]) * 100
        rows[m] = (len(p), q)
        print(f"   {m:>9,} {len(p):>7,} {q.iloc[0]:>7.2f}% {q.iloc[1]:>7.2f}% {q.iloc[2]:>7.2f}%")
    # largest cutoff leaving >= 2,000 genes (NTR falls with expression, so a very high
    # cutoff keeps only stable genes); fall back to >= 300 genes
    cand = [m for m, (n, _) in rows.items() if n >= 2000] or \
           [m for m, (n, _) in rows.items() if n >= 300]
    min_umis = max(cand) if cand else (min(rows) if rows else 50)
    pool = df[(df.UMIs >= min_umis) & (df.cells >= min_cells)]
    q10, q50, q90 = (pool.NTR.quantile([0.1, 0.5, 0.9]) * 100).values

    hi = pool[pool.NTR * 100 >= q90]
    lq = hi.labeled.quantile([0.25, 0.5]).values
    print(f"\n3. High-NTR genes (top 10% at >= {min_umis} UMIs, {len(hi)} genes): labeled UMIs "
          f"25th pct {lq[0]:.0f}, median {lq[1]:.0f}")
    min_lab = int(max(10, min(50, round(lq[0], -1))))

    top = pool.sort_values("UMIs", ascending=False).head(200).index
    per = []
    for gi in top:
        col = np.asarray(tot[:, gi].todense()).ravel()
        col = col[col > 0]
        if len(col):
            per.append(np.median(col))
    per = np.array(per) if per else np.array([1])
    print(f"\n4. UMIs per expressing cell, median over the 200 best-expressed genes: "
          f"{np.median(per):.1f} (25th pct {np.percentile(per, 25):.1f}, "
          f"75th pct {np.percentile(per, 75):.1f})")
    # high-NTR genes are often less expressed: keep the cell filter modest
    min_cell_umis = int(max(3, min(5, round(np.percentile(per, 25)))))

    print("\n5. Suggested defaults")
    sug = [
        ("--min-gene-umis", min_umis, "largest cutoff that still leaves >= 2,000 genes"),
        ("--min-cells", min_cells, "max(20, 1% of cells)"),
        ("--min-gene-labeled", min_lab, "about the 25th pct of labeled UMIs among high-NTR genes"),
        ("--high-ntr-percentile", 90, f"top 10% (NTR >= {q90:.2f}% here)"),
        ("--low-ntr-percentile", 10, f"bottom 10% (NTR <= {q10:.2f}% here)"),
        ("--rank-by", "umis", "best-expressed genes among those passing -> fuller panels"),
        ("--min-cell-umis", min_cell_umis, "modest (<= 5): high-NTR genes are less expressed"),
        ("--min-cell-labeled", 2, "at least 2 labeled transcripts of the high-NTR gene"),
        ("--min-cell-ntr", round(2 / 3 * q90, 1), "2/3 of the gene-level 90th pct (cell NTR is noisy)"),
        ("--max-cell-ntr", round(max(q10, 0.0), 1), "low-NTR gene in the cell at most the 10th pct"),
    ]
    for opt, val, why in sug:
        print(f"   {opt:22s} {str(val):>8}   {why}")
    if overall < 0.02:
        print("\n   NOTE: overall NTR is very low (unlabeled or weakly labeled sample). Use the "
              "percentile options; fixed NTR cutoffs will leave few genes.")
    opts = " ".join(f"{o} {v}" for o, v, _ in sug)
    print("\nSuggested command:")
    print(f"python scripts/plot_tc_panel.py --bam \"$bam\" --fasta 00_input/genome.fa "
          f"--gtf 00_input/annotation.gtf --adata {a.adata} --no-mt-ribo \\\n"
          f"    --high-ntr-genes 1 --low-ntr-genes 1 {opts} \\\n"
          f"    --snp-csv 01_mask/DONOR_snps_0based.csv --vmax 4 "
          f"--out 04_plots/high_vs_low_NTR.pdf")


if __name__ == "__main__":
    main()
