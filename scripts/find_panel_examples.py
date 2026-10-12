#!/usr/bin/env python3
"""
find_panel_examples.py -- search your data for the best example for the
scNT-seq Fig. 1c-style panel: ONE cell in which a fast-turnover gene has many
labeled transcripts and a stable gene has many mostly unlabeled transcripts.

For every candidate (fast gene, stable gene, cell) it uses that cell's UMI
counts from the dynast count output (labeled_TC / unlabeled_TC layers) and
scores
    score = (labeled fraction fast - labeled fraction stable) * sqrt(min UMIs)
so a good example has a strong contrast AND enough transcripts in both panels.
It prints the top examples with a ready-to-run plot_tc_panel.py command.

Usage:
  python find_panel_examples.py \
      --adata 02_dynast_count/TC_rm_bg_0based/adata.h5ad \
      --bam /path/split_bams/<donor>.bam --snp-csv <donor>_snps_0based.csv \
      --fasta 00_input/genome.fa --gtf 00_input/annotation.gtf
  # force genes:  --fast-genes FOS JUN  --stable-genes EEF1A1 RPLP0
"""
import argparse
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp


def get_layers(ad, conv):
    for lab, unl in ((f"labeled_{conv}", f"unlabeled_{conv}"), (f"X_l_{conv}", f"X_n_{conv}")):
        if lab in ad.layers and unl in ad.layers:
            return sp.csc_matrix(ad.layers[lab]), sp.csc_matrix(ad.layers[unl])
    sys.exit(f"ERROR: no labeled/unlabeled {conv} layers; found: {list(ad.layers.keys())}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--adata", required=True, help="dynast count adata.h5ad (masked 0-based run)")
    ap.add_argument("--conversion", default="TC")
    ap.add_argument("--fast-genes", nargs="*", default=[], help="fast candidates (default: auto)")
    ap.add_argument("--stable-genes", nargs="*", default=[], help="stable candidates (default: auto)")
    ap.add_argument("--n-candidates", type=int, default=15, help="auto candidates per group")
    ap.add_argument("--min-cells", type=int, default=30, help="auto: gene expressed in >= N cells")
    ap.add_argument("--top-expressed", type=int, default=1000, help="auto: pick from N most expressed")
    ap.add_argument("--no-mt-ribo", action="store_true", help="auto: skip MT-, RPL*, RPS* genes")
    ap.add_argument("--min-umis", type=int, default=15,
                    help="each gene needs >= N UMIs in the cell (panel readability)")
    ap.add_argument("--max-umis", type=int, default=80,
                    help="prefer cells where a gene has <= N UMIs (larger panels get crowded)")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--out", default="panel_examples.tsv")
    # only used to print ready-to-run commands
    ap.add_argument("--bam", default="<pre-consensus BAM>")
    ap.add_argument("--fasta", default="00_input/genome.fa")
    ap.add_argument("--gtf", default="00_input/annotation.gtf")
    ap.add_argument("--snp-csv", default="")
    a = ap.parse_args()

    import anndata
    ad = anndata.read_h5ad(a.adata)
    lab, unl = get_layers(ad, a.conversion)
    tot = lab + unl
    names = (ad.var["gene_name"].astype(str).values if "gene_name" in ad.var
             else ad.var_names.astype(str).values)
    cells = ad.obs_names.astype(str).values

    g_tot = np.asarray(tot.sum(axis=0)).ravel()
    g_lab = np.asarray(lab.sum(axis=0)).ravel()
    g_cells = np.asarray((tot > 0).sum(axis=0)).ravel()
    g_frac = np.divide(g_lab, g_tot, out=np.full_like(g_tot, np.nan, dtype=float), where=g_tot > 0)
    overall = g_lab.sum() / g_tot.sum()
    print(f"{ad.n_obs:,} cells x {ad.n_vars:,} genes; overall labeled fraction {overall:.1%}")
    if overall < 0.02:
        print("NOTE: very low labeling -- if this sample had no 4sU, there is no good example "
              "here. Use a labeled sample.")

    order = np.argsort(-g_tot)
    pool = [i for i in order[: a.top_expressed] if g_cells[i] >= a.min_cells]
    if a.no_mt_ribo:
        pool = [i for i in pool if not pd.Series(names[i]).str.match(r"^(MT-|RPL|RPS)").iloc[0]]

    def pick(given, fast):
        if given:
            idx = []
            for g in given:
                w = np.where(names == g)[0]
                if len(w) == 0:
                    print(f"WARNING: {g} not in adata")
                else:
                    idx.append(w[0])
            return idx
        srt = sorted(pool, key=lambda i: g_frac[i], reverse=fast)
        return srt[: a.n_candidates]

    fast_idx = pick(a.fast_genes, True)
    stab_idx = pick(a.stable_genes, False)
    print("fast candidates  : " + ", ".join(f"{names[i]} ({g_frac[i]:.0%})" for i in fast_idx))
    print("stable candidates: " + ", ".join(f"{names[i]} ({g_frac[i]:.0%})" for i in stab_idx))

    def dense(m, i):
        return np.asarray(m[:, i].todense()).ravel()

    rows = []
    for fi in fast_idx:
        nf, lf = dense(tot, fi), dense(lab, fi)
        for si in stab_idx:
            if si == fi:
                continue
            ns, ls = dense(tot, si), dense(lab, si)
            ok = (nf >= a.min_umis) & (ns >= a.min_umis)
            for c in np.where(ok)[0]:
                ff, fs = lf[c] / nf[c], ls[c] / ns[c]
                m = min(nf[c], ns[c])
                crowd = max(1.0, max(nf[c], ns[c]) / a.max_umis)   # penalise huge panels
                rows.append((names[fi], names[si], cells[c], int(nf[c]), int(lf[c]), ff,
                             int(ns[c]), int(ls[c]), fs, (ff - fs) * np.sqrt(m) / crowd))
    if not rows:
        sys.exit(f"No cell has >= {a.min_umis} UMIs for a fast AND a stable candidate. "
                 "Lower --min-umis, or give other genes / another sample.")

    df = pd.DataFrame(rows, columns=["fast_gene", "stable_gene", "cell", "fast_UMIs",
                                     "fast_labeled", "fast_frac", "stable_UMIs",
                                     "stable_labeled", "stable_frac", "score"])
    df = df.sort_values("score", ascending=False)
    df.to_csv(a.out, sep="\t", index=False, float_format="%.4g")

    best_pairs = df.drop_duplicates(["fast_gene", "stable_gene"]).head(a.top)
    fmt = {"fast_frac": "{:.0%}".format, "stable_frac": "{:.0%}".format, "score": "{:.2f}".format}
    print(f"\n=== Best examples (best cell per gene pair; all combinations in {a.out})")
    print(best_pairs.to_string(index=False, formatters=fmt))

    print("\n=== Commands for the top 3")
    for _, r in best_pairs.head(3).iterrows():
        snp = f" --snp-csv {a.snp_csv}" if a.snp_csv else ""
        print(f"python scripts/plot_tc_panel.py --bam {a.bam} --fasta {a.fasta} --gtf {a.gtf} "
              f"--genes {r.fast_gene} {r.stable_gene} --cell {r.cell}{snp} "
              f"--max-umis {max(r.fast_UMIs, r.stable_UMIs)} "
              f"--out 04_plots/panel_{r.fast_gene}_{r.stable_gene}_{r.cell}.pdf")


if __name__ == "__main__":
    main()
