#!/usr/bin/env python3
"""
rank_turnover_genes.py -- pick genes for the scNT-seq Fig. 1c-style panel from
YOUR data: well-expressed stable genes (low labeled fraction, "housekeeping"
panel) and fast-turnover genes (high labeled fraction), using the dynast count
output of the corrected (0-based) masked run.

Labeled fraction of a gene = labeled UMIs / all UMIs (summed over cells), from
the dynast layers labeled_TC / unlabeled_TC (or X_l_TC / X_n_TC).

With --pair-with FOS, it also counts how many cells express each candidate
TOGETHER with FOS, so plot_tc_panel.py can find one cell for both panels.

Usage:
  python rank_turnover_genes.py \
      --adata 02_dynast_count/TC_rm_bg_0based/adata.h5ad \
      --pair-with FOS --top 20 --out 03_checks/gene_turnover_ranking.tsv
"""
import argparse
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

KNOWN_HOUSEKEEPING = ["EEF1A1", "RPLP0", "RPL13A", "RPS18", "GAPDH", "ACTB", "PPIA",
                      "B2M", "TUBB", "UBC", "HPRT1", "TBP", "PGK1", "YWHAZ", "SDHA", "MAPT"]


def get_layers(ad, conv):
    for lab, unl in ((f"labeled_{conv}", f"unlabeled_{conv}"), (f"X_l_{conv}", f"X_n_{conv}")):
        if lab in ad.layers and unl in ad.layers:
            return ad.layers[lab], ad.layers[unl], (lab, unl)
    sys.exit(f"ERROR: no labeled/unlabeled {conv} layers; found: {list(ad.layers.keys())}")


def colsum(m):
    return np.asarray(m.sum(axis=0)).ravel()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--adata", required=True, help="dynast count adata.h5ad (masked 0-based run)")
    ap.add_argument("--conversion", default="TC")
    ap.add_argument("--pair-with", nargs="*", default=[],
                    help="gene(s) that must be in the same cell (e.g. FOS)")
    ap.add_argument("--min-cells", type=int, default=50,
                    help="gene must have >=1 UMI in at least this many cells")
    ap.add_argument("--top-expressed", type=int, default=500,
                    help="stable/fast lists are drawn from the N most expressed genes")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--no-mt-ribo", action="store_true",
                    help="exclude MT- and RPL/RPS genes from the lists")
    ap.add_argument("--out", default="gene_turnover_ranking.tsv")
    a = ap.parse_args()

    import anndata
    ad = anndata.read_h5ad(a.adata)
    lab, unl, names = get_layers(ad, a.conversion)
    lab, unl = sp.csr_matrix(lab), sp.csr_matrix(unl)
    tot = lab + unl
    names_col = ad.var["gene_name"].astype(str).values if "gene_name" in ad.var else ad.var_names.astype(str).values
    print(f"{ad.n_obs:,} cells x {ad.n_vars:,} genes; layers: {names[0]}, {names[1]}")

    df = pd.DataFrame({
        "gene": names_col,
        "gene_id": ad.var_names.astype(str),
        "total_UMIs": colsum(tot).astype(np.int64),
        "labeled_UMIs": colsum(lab).astype(np.int64),
        "cells_expressing": np.asarray((tot > 0).sum(axis=0)).ravel(),
    })
    df["labeled_fraction"] = np.where(df.total_UMIs > 0, df.labeled_UMIs / df.total_UMIs, np.nan)
    df["UMIs_per_expressing_cell"] = np.where(df.cells_expressing > 0,
                                              df.total_UMIs / df.cells_expressing, 0)

    # co-expression with the partner gene(s)
    partner_cells = None
    for g in a.pair_with:
        idx = np.where(names_col == g)[0]
        if len(idx) == 0:
            print(f"WARNING: --pair-with {g} not found in adata")
            continue
        cells_g = np.asarray(tot[:, idx].sum(axis=1)).ravel() > 0
        print(f"{g}: {cells_g.sum():,} cells with >=1 UMI, {int(tot[:, idx].sum()):,} UMIs, "
              f"labeled fraction {df.loc[idx[0], 'labeled_fraction']:.1%}")
        partner_cells = cells_g if partner_cells is None else (partner_cells & cells_g)
    if partner_cells is not None:
        df["cells_with_partner"] = np.asarray((tot[partner_cells] > 0).sum(axis=0)).ravel()

    df = df.sort_values("total_UMIs", ascending=False).reset_index(drop=True)
    df["expression_rank"] = np.arange(1, len(df) + 1)
    df.to_csv(a.out, sep="\t", index=False, float_format="%.4g")

    overall = df.labeled_UMIs.sum() / df.total_UMIs.sum()
    print(f"Overall labeled fraction (all genes): {overall:.1%}")
    if overall < 0.02:
        print("NOTE: very low overall labeling -- if this sample had no 4sU, labeled fractions are "
              "background only and cannot rank turnover. Use a labeled sample.")

    pool = df[(df.cells_expressing >= a.min_cells)].head(a.top_expressed)
    if a.no_mt_ribo:
        pool = pool[~pool.gene.str.match(r"^(MT-|RPL|RPS)")]
    cols = ["gene", "expression_rank", "total_UMIs", "cells_expressing",
            "UMIs_per_expressing_cell", "labeled_fraction"]
    if "cells_with_partner" in df:
        cols.append("cells_with_partner")
    fmt = {"labeled_fraction": "{:.1%}".format, "UMIs_per_expressing_cell": "{:.1f}".format}

    print(f"\n=== STABLE candidates (lowest labeled fraction among the {a.top_expressed} most "
          f"expressed genes in >= {a.min_cells} cells) -> bottom panel")
    print(pool.nsmallest(a.top, "labeled_fraction")[cols].to_string(index=False, formatters=fmt))
    print(f"\n=== FAST candidates (highest labeled fraction) -> top panel")
    print(pool.nlargest(a.top, "labeled_fraction")[cols].to_string(index=False, formatters=fmt))

    known = df[df.gene.isin(KNOWN_HOUSEKEEPING)]
    if len(known):
        print("\n=== Common housekeeping genes in this sample")
        print(known[cols].to_string(index=False, formatters=fmt))
    print(f"\nFull table: {a.out}")


if __name__ == "__main__":
    main()
