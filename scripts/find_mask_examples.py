#!/usr/bin/env python3
"""
find_mask_examples.py -- find the best (gene, cell) examples for a before/after
SNP-mask plot that show BOTH effects in one cell:
  * many molecules whose "label" came only from the SNP (removed by the mask)
  * many molecules that stay labeled after masking (real 4sU T>C)

Uses the dynast count output without and with the 0-based mask (same UMIs, only
the labeled/unlabeled split differs). Per gene and cell:
    removed  = labeled UMIs without mask - labeled UMIs with mask
    kept     = labeled UMIs with mask   (real T>C remaining)
    score    = sqrt(removed * kept), penalised when the panel gets very large

Prints the best cell per gene and ready-to-run plot_tc_panel.py commands.

Usage:
  python find_mask_examples.py \
      --unmasked 02_dynast_count/TC/adata.h5ad \
      --masked   02_dynast_count/TC_rm_bg_0based/adata.h5ad \
      --snp-csv  01_mask/<donor>_snps_0based.csv
"""
import argparse
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp


def load(path, conv):
    import anndata
    ad = anndata.read_h5ad(path)
    for lab, unl in ((f"labeled_{conv}", f"unlabeled_{conv}"), (f"X_l_{conv}", f"X_n_{conv}")):
        if lab in ad.layers and unl in ad.layers:
            names = (ad.var["gene_name"].astype(str).values if "gene_name" in ad.var
                     else ad.var_names.astype(str).values)
            return ad, sp.csr_matrix(ad.layers[lab]), sp.csr_matrix(ad.layers[unl]), names
    sys.exit(f"ERROR: no labeled/unlabeled {conv} layers in {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--unmasked", required=True, help="adata.h5ad WITHOUT mask")
    ap.add_argument("--masked", required=True, help="adata.h5ad WITH the 0-based mask")
    ap.add_argument("--conversion", default="TC")
    ap.add_argument("--min-removed", type=int, default=2,
                    help="cell needs >= N molecules whose label was removed by the mask")
    ap.add_argument("--min-kept", type=int, default=3,
                    help="cell needs >= N molecules still labeled after masking")
    ap.add_argument("--max-umis", type=int, default=80,
                    help="prefer cells with <= N UMIs of the gene (readable panels)")
    ap.add_argument("--no-mt-ribo", action="store_true", help="skip MT-, RPL*, RPS* genes")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--out", default="mask_examples.tsv")
    ap.add_argument("--bam", default="00_input/consensus.bam")
    ap.add_argument("--fasta", default="00_input/genome.fa")
    ap.add_argument("--gtf", default="00_input/annotation.gtf")
    ap.add_argument("--snp-csv", default="01_mask/<donor>_snps_0based.csv")
    ap.add_argument("--plot-args", default="--min-bq 28 --control-shift 1 --window 150")
    a = ap.parse_args()

    ad0, lab0, unl0, names = load(a.unmasked, a.conversion)
    ad1, lab1, unl1, _ = load(a.masked, a.conversion)
    cells = ad0.obs_names.intersection(ad1.obs_names)
    genes = ad0.var_names.intersection(ad1.var_names)
    r0, r1 = ad0.obs_names.get_indexer(cells), ad1.obs_names.get_indexer(cells)
    g0, g1 = ad0.var_names.get_indexer(genes), ad1.var_names.get_indexer(genes)
    l0 = lab0[r0][:, g0]
    tot = (lab0 + unl0)[r0][:, g0]
    l1 = lab1[r1][:, g1]
    gname = names[g0]

    overall0 = l0.sum() / max(tot.sum(), 1)
    overall1 = l1.sum() / max(tot.sum(), 1)
    print(f"{len(cells):,} cells, {len(genes):,} genes; labeled fraction "
          f"{overall0:.2%} without mask -> {overall1:.2%} with mask")
    if overall1 < 0.02:
        print("NOTE: very little labeling remains after masking. This sample probably had little "
              "or no 4sU, so few cells will have many real T>C. Use a labeled sample (later day).")

    rem = sp.coo_matrix(l0 - l1)
    keep = rem.data >= a.min_removed
    ci, gi, rv = rem.row[keep], rem.col[keep], rem.data[keep]
    if len(ci) == 0:
        sys.exit(f"No cell has >= {a.min_removed} molecules whose label was removed. "
                 "Lower --min-removed.")
    l1d = np.asarray(l1[ci, gi]).ravel()
    tod = np.asarray(tot[ci, gi]).ravel()
    df = pd.DataFrame({"gene": gname[gi], "cell": np.asarray(cells)[ci],
                       "UMIs": tod.astype(int), "labeled_no_mask": (l1d + rv).astype(int),
                       "removed_by_mask": rv.astype(int), "labeled_with_mask": l1d.astype(int)})
    df = df[df.labeled_with_mask >= a.min_kept]
    if a.no_mt_ribo:
        df = df[~df.gene.str.match(r"^(MT-|RPL|RPS)")]
    if df.empty:
        sys.exit(f"No cell has >= {a.min_removed} removed AND >= {a.min_kept} remaining labeled "
                 "molecules in one gene. Lower --min-kept / --min-removed, or use a labeled sample.")
    crowd = np.maximum(1.0, df.UMIs / a.max_umis)
    df["score"] = np.sqrt(df.removed_by_mask * df.labeled_with_mask) / crowd
    df["frac_with_mask"] = df.labeled_with_mask / df.UMIs
    df = df.sort_values("score", ascending=False)
    df.to_csv(a.out, sep="\t", index=False, float_format="%.4g")

    best = df.drop_duplicates("gene").head(a.top)
    fmt = {"score": "{:.2f}".format, "frac_with_mask": "{:.0%}".format}
    print(f"\n=== Best cell per gene (all candidates in {a.out})")
    print(best.to_string(index=False, formatters=fmt))

    print("\n=== Commands for the top 5")
    for _, r in best.head(5).iterrows():
        print(f"python scripts/plot_tc_panel.py --bam {a.bam} --fasta {a.fasta} --gtf {a.gtf} "
              f"--genes {r.gene} --before-after --cell {r.cell} --snp-csv {a.snp_csv} "
              f"{a.plot_args} --out 04_plots/before_after_{r.gene}_{r.cell}.pdf")


if __name__ == "__main__":
    main()
