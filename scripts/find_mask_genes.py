#!/usr/bin/env python3
"""
find_mask_genes.py -- find candidate genes for the before/after SNP-mask plot:
genes whose labeled (T>C) UMIs drop the most when the WGS SNP mask is applied.

Compares the dynast count output WITHOUT a mask (TC/adata.h5ad) and WITH the
0-based mask (TC_rm_bg_0based/adata.h5ad). Both runs count the same UMIs; only
the labeled/unlabeled split changes, so for each gene:
    removed = labeled UMIs without mask - labeled UMIs with mask
and the cell where most labeled UMIs were removed is reported.

Optionally joins the site table from find_snp_examples.py (--sites) to give the
SNP position for each gene, and prints ready-to-run plot_tc_panel.py commands.

Usage:
  python find_mask_genes.py \
      --unmasked 02_dynast_count/TC/adata.h5ad \
      --masked   02_dynast_count/TC_rm_bg_0based/adata.h5ad \
      --sites    03_checks/snp_examples_DONOR.tsv \
      --snp-csv 01_mask/<donor>_snps_0based.csv      (--bam defaults to 00_input/consensus.bam)
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
            return ad, sp.csc_matrix(ad.layers[lab]), sp.csc_matrix(ad.layers[unl]), names
    sys.exit(f"ERROR: no labeled/unlabeled {conv} layers in {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--unmasked", required=True, help="adata.h5ad of the run WITHOUT mask")
    ap.add_argument("--masked", required=True, help="adata.h5ad of the run WITH the 0-based mask")
    ap.add_argument("--conversion", default="TC")
    ap.add_argument("--sites", help="snp_examples TSV from find_snp_examples.py (optional)")
    ap.add_argument("--genes", nargs="*", default=[], help="only report these genes (e.g. FOS JUN)")
    ap.add_argument("--min-umis", type=int, default=50, help="gene needs >= N UMIs in total")
    ap.add_argument("--no-mt-ribo", action="store_true", help="skip MT-, RPL*, RPS* genes")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--out", default="mask_affected_genes.tsv")
    ap.add_argument("--bam", default="00_input/consensus.bam",
                    help="BAM for the printed commands (consensus = what dynast counted)")
    ap.add_argument("--fasta", default="00_input/genome.fa")
    ap.add_argument("--gtf", default="00_input/annotation.gtf")
    ap.add_argument("--snp-csv", default="01_mask/<donor>_snps_0based.csv")
    ap.add_argument("--plot-args", default="--min-bq 28 --control-shift 1 --window 50",
                    help="extra plot_tc_panel.py options in the printed commands")
    a = ap.parse_args()

    ad0, lab0, unl0, names0 = load(a.unmasked, a.conversion)
    ad1, lab1, unl1, names1 = load(a.masked, a.conversion)

    # align cells and genes of the two runs
    cells = ad0.obs_names.intersection(ad1.obs_names)
    genes = ad0.var_names.intersection(ad1.var_names)
    i0c, i1c = ad0.obs_names.get_indexer(cells), ad1.obs_names.get_indexer(cells)
    i0g, i1g = ad0.var_names.get_indexer(genes), ad1.var_names.get_indexer(genes)
    lab0, tot0 = lab0[i0c][:, i0g], (lab0 + unl0)[i0c][:, i0g]
    lab1 = lab1[i1c][:, i1g]
    gname = names0[i0g]
    print(f"{len(cells):,} cells and {len(genes):,} genes in both runs")

    removed = sp.csc_matrix(lab0 - lab1)
    neg = removed.min() if removed.nnz else 0
    if neg < 0:
        print("NOTE: some cells/genes have MORE labeled UMIs with the mask "
              "(e.g. the masked run also used --snp-threshold, or different settings)")

    tot = np.asarray(tot0.sum(axis=0)).ravel()
    l0 = np.asarray(lab0.sum(axis=0)).ravel()
    l1 = np.asarray(lab1.sum(axis=0)).ravel()
    rem = l0 - l1
    cells_aff = np.asarray((removed > 0).sum(axis=0)).ravel()
    best_idx = np.asarray(removed.argmax(axis=0)).ravel()
    best_n = np.asarray(removed.max(axis=0).todense()).ravel()

    df = pd.DataFrame({
        "gene": gname, "gene_id": genes.astype(str), "total_UMIs": tot.astype(int),
        "labeled_no_mask": l0.astype(int), "labeled_with_mask": l1.astype(int),
        "removed_labeled": rem.astype(int),
        "frac_no_mask": np.divide(l0, tot, out=np.zeros_like(l0, float), where=tot > 0),
        "frac_with_mask": np.divide(l1, tot, out=np.zeros_like(l1, float), where=tot > 0),
        "cells_affected": cells_aff.astype(int),
        "best_cell": np.asarray(cells)[best_idx], "best_cell_removed": best_n.astype(int),
    })
    df["frac_drop"] = df.frac_no_mask - df.frac_with_mask
    df.loc[df.best_cell_removed <= 0, "best_cell"] = ""

    if a.sites:
        st = pd.read_csv(a.sites, sep="\t", comment=None)
        st = st[pd.to_numeric(st["rank"], errors="coerce").notna()]
        if "counted_as_TC" in st:
            st = st[st.counted_as_TC == "yes"]
        top_site = st.sort_values("events", ascending=False).drop_duplicates("gene")
        df = df.merge(top_site[["gene", "chrom", "pos_1based", "events"]]
                      .rename(columns={"events": "site_TC_events"}), on="gene", how="left")

    df = df.sort_values(["removed_labeled", "frac_drop"], ascending=False)
    df.to_csv(a.out, sep="\t", index=False, float_format="%.4g")

    view = df[df.total_UMIs >= a.min_umis]
    if a.no_mt_ribo:
        view = view[~view.gene.str.match(r"^(MT-|RPL|RPS)")]
    if a.genes:
        view = df[df.gene.isin(a.genes)]
        missing = set(a.genes) - set(view.gene)
        if missing:
            print(f"not found: {', '.join(sorted(missing))}")
    else:
        view = view[view.removed_labeled > 0]

    cols = ["gene", "total_UMIs", "labeled_no_mask", "labeled_with_mask", "removed_labeled",
            "frac_no_mask", "frac_with_mask", "cells_affected", "best_cell", "best_cell_removed"]
    if "pos_1based" in df:
        cols += ["pos_1based"]
    fmt = {"frac_no_mask": "{:.1%}".format, "frac_with_mask": "{:.1%}".format,
           "pos_1based": lambda v: "" if pd.isna(v) else f"{int(v):,}"}
    print(f"\n=== Genes whose labeled UMIs change most with the mask (total UMIs >= {a.min_umis})")
    print(view.head(a.top)[cols].to_string(index=False, formatters=fmt, na_rep="") if len(view) else "  none")
    print(f"\nTotal labeled UMIs removed by the mask: {int(rem.sum()):,} of {int(l0.sum()):,} "
          f"({rem.sum() / max(l0.sum(), 1):.1%}); genes affected: {(rem > 0).sum():,}")
    print(f"Full table: {a.out}")

    print("\n=== Commands for the top 5 (before/after plot; site and window are chosen automatically)")
    for _, r in view[view.removed_labeled > 0].head(5).iterrows():
        print(f"python scripts/plot_tc_panel.py --bam {a.bam} --fasta {a.fasta} --gtf {a.gtf} "
              f"--genes {r.gene} --before-after --snp-csv {a.snp_csv} {a.plot_args} "
              f"--out 04_plots/before_after_{r.gene}_consensus.pdf")


if __name__ == "__main__":
    main()
