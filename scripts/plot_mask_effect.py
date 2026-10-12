#!/usr/bin/env python3
"""
plot_mask_effect.py -- compare per-cell T>C rates between dynast count runs
(e.g. no mask, old 1-based mask, new 0-based mask).

Reads <run_dir>/rates.csv (TC column = per-cell T>C rate, strand-corrected).
The FIRST run is the reference (normally the unmasked run).

Usage:
  python plot_mask_effect.py \
      --run "No mask=donor_count/<donor>/TC" \
      --run "Old mask (1-based)=donor_count/<donor>/TC_rm_bg" \
      --run "New mask (0-based)=donor_count/<donor>/TC_rm_bg_0based" \
      --title "<donor>" --out mask_effect_<donor>.pdf
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_tc(run_dir):
    path = os.path.join(run_dir, "rates.csv")
    if not os.path.exists(path):
        sys.exit(f"ERROR: {path} not found (did dynast count finish? check its log)")
    df = pd.read_csv(path)
    return df.set_index("barcode")["TC"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True,
                    help='"label=run_dir"; give 2 or more, the first is the reference')
    ap.add_argument("--title", default="")
    ap.add_argument("--out", default="mask_effect.pdf")
    a = ap.parse_args()

    labels, series = [], []
    for item in a.run:
        label, d = item.split("=", 1)
        labels.append(label)
        series.append(read_tc(d))
    df = pd.concat(series, axis=1, keys=labels, join="inner").dropna()
    if df.empty:
        sys.exit("ERROR: no barcodes shared between runs")

    ref = labels[0]
    print(f"{len(df):,} cells in common")
    print(f"{'run':30s} {'median TC rate':>15s} {'vs ' + ref:>15s} {'cells lower':>12s}")
    for lab in labels:
        m = df[lab].median()
        change = m / df[ref].median() - 1 if df[ref].median() > 0 else np.nan
        lower = (df[lab] < df[ref]).mean()
        print(f"{lab:30s} {m:15.6f} {change:+15.1%} {lower:12.1%}")

    n_other = len(labels) - 1
    fig, axes = plt.subplots(1, 1 + n_other, figsize=(4.2 * (1 + n_other), 4))
    axes = np.atleast_1d(axes)

    # (1) distribution of per-cell T>C rate in each run
    ax = axes[0]
    ax.boxplot([df[l].values for l in labels], showfliers=False, widths=0.6)
    for i, l in enumerate(labels, 1):
        x = np.random.default_rng(0).normal(i, 0.06, len(df))
        ax.scatter(x, df[l].values, s=3, alpha=0.3, color="0.3")
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Per-cell T>C rate")
    ax.set_title("Per-cell T>C rate")

    # (2) per-cell scatter vs the reference run
    lim = [0, float(df.values.max()) * 1.05]
    for ax, lab in zip(axes[1:], labels[1:]):
        ax.scatter(df[ref], df[lab], s=4, alpha=0.4)
        ax.plot(lim, lim, ls="--", color="red", lw=1, label="y = x (no change)")
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel(f"T>C rate: {ref}")
        ax.set_ylabel(f"T>C rate: {lab}")
        drop = 1 - df[lab].median() / df[ref].median()
        ax.set_title(f"{lab}\nmedian {drop:.1%} lower")
        ax.legend(frameon=False, fontsize=8, loc="upper left")

    fig.suptitle(a.title)
    fig.tight_layout()
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    print(f"\nSaved: {a.out}")


if __name__ == "__main__":
    main()
