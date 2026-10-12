#!/usr/bin/env python3
"""
check_snp_mask.py -- verify that a WGS SNP mask is correct for dynast and that
dynast actually applied it.

Checks (each runs only if its inputs are given):
  1. FORMAT     header, columns, conversion labels, contig names     (--mask [--bam])
  2. COORDINATE reference base at genome_i must be T (TC) / A (AG)   (--fasta)
                -> proves the file is 0-based like dynast
  3. DYNAST     conversions in an UNMASKED run land exactly on mask
                positions, not 1 bp away                             (--conversions)
  4. EFFECT     per-cell T>C rate drops in the masked run            (--unmasked --masked)

Usage:
  python check_snp_mask.py --mask snps_union_all.csv --fasta GRCh38.fa \
      --bam consensus.bam --conversions TC/conversions.csv \
      --unmasked TC --masked TC_rm_bg
"""
import argparse
import os
import random
import sys
from collections import defaultdict

PASS, FAIL, WARN = "PASS", "FAIL", "WARN"
results = []


def report(name, status, msg):
    results.append((name, status))
    print(f"[{status}] {name}: {msg}")


# ---------------------------------------------------------------------------
def load_mask(path):
    """Return header, {conversion: {contig: set(pos)}}, n_rows, problems."""
    mask = defaultdict(lambda: defaultdict(set))
    problems, n = [], 0
    with open(path) as f:
        header = f.readline().strip()
        for line in f:
            p = line.rstrip("\n").split(",")
            n += 1
            if len(p) != 3 or not p[1].isdigit():
                if len(problems) < 5:
                    problems.append(line.strip())
                continue
            mask[p[2]][p[0]].add(int(p[1]))
    return header, mask, n, problems


def check_format(header, mask, n, problems, bam):
    if header != "contig,genome_i,conversion":
        report("FORMAT header", FAIL,
               f"first line is '{header}'; dynast needs exactly 'contig,genome_i,conversion' "
               "(otherwise a 3-column file raises 'Failed to parse')")
    else:
        report("FORMAT header", PASS, "contig,genome_i,conversion")

    if problems:
        report("FORMAT rows", FAIL, f"malformed rows, e.g. {problems}")
    else:
        report("FORMAT rows", PASS, f"{n:,} rows, all 3 columns with integer genome_i")

    labels = set(mask)
    if labels - {"TC", "AG"}:
        report("FORMAT conversion", WARN, f"unexpected labels {sorted(labels)} (expected TC, AG)")
    else:
        counts = {c: sum(len(s) for s in mask[c].values()) for c in labels}
        report("FORMAT conversion", PASS, f"site counts {counts}")

    if bam:
        import pysam
        refs = set(pysam.AlignmentFile(bam, "rb").references)
        contigs = {c for d in mask.values() for c in d}
        missing = sorted(contigs - refs)
        if missing == sorted(contigs):
            report("FORMAT contigs", FAIL, f"no mask contig is in the BAM header (mask: {sorted(contigs)[:5]}, BAM: {sorted(refs)[:5]})")
        elif missing:
            report("FORMAT contigs", WARN, f"contigs not in BAM: {missing[:10]}")
        else:
            report("FORMAT contigs", PASS, "all mask contigs are present in the BAM")


# ---------------------------------------------------------------------------
def check_coordinates(mask, fasta, n_sample, seed=1):
    """At a 0-based position the reference base must equal the SNP REF base."""
    import pysam
    fa = pysam.FastaFile(fasta)
    expect = {"TC": "T", "AG": "A"}
    rows = [(conv, c, p) for conv, d in mask.items() if conv in expect
            for c, s in d.items() if c in fa.references for p in s]
    if not rows:
        report("COORDINATE", FAIL, "no mask rows on contigs present in the FASTA")
        return
    random.Random(seed).shuffle(rows)
    rows = rows[:n_sample]
    hit0 = hit1 = 0
    for conv, c, p in rows:
        b = fa.fetch(c, max(p - 1, 0), p + 2).upper()   # bases at p-1, p, p+1 (0-based)
        if p == 0:
            b = "N" + b
        hit0 += b[1] == expect[conv]                      # file is 0-based
        hit1 += b[0] == expect[conv]                      # file is 1-based (VCF POS)
    f0, f1 = hit0 / len(rows), hit1 / len(rows)
    msg = (f"ref base matches at genome_i (0-based): {f0:.1%}; "
           f"at genome_i-1 (i.e. file is 1-based): {f1:.1%}  [n={len(rows):,}]")
    if f0 >= 0.99:
        report("COORDINATE", PASS, msg + " -> file is 0-based, matches dynast")
    elif f1 >= 0.99:
        report("COORDINATE", FAIL, msg + " -> file holds VCF POS (1-based); subtract 1")
    else:
        report("COORDINATE", FAIL, msg + " -> wrong reference build, strand or columns?")


# ---------------------------------------------------------------------------
def check_dynast(mask, conversions, quality, max_lines):
    """Count TC/AG conversions in an UNMASKED dynast run at mask positions vs +/-1."""
    if not os.path.exists(conversions):
        report("DYNAST", FAIL, f"{conversions} not found -- did dynast count finish?")
        return
    exact = plus1 = minus1 = total = 0
    with open(conversions) as f:
        first = f.readline()
        if not first.startswith("read_id"):
            f.seek(0)
        for i, line in enumerate(f):
            if max_lines and i >= max_lines:
                break
            p = line.rstrip("\n").split(",")
            conv = p[4]
            if conv not in ("TC", "AG") or int(p[5]) <= quality:
                continue
            total += 1
            s = mask.get(conv, {}).get(p[2])
            if not s:
                continue
            g = int(p[3])
            exact += g in s
            plus1 += (g + 1) in s
            minus1 += (g - 1) in s
    if total == 0:
        size = os.path.getsize(conversions)
        report("DYNAST", FAIL, f"no TC/AG conversions with quality > {quality} in {conversions} "
               f"({size:,} bytes) -- did dynast count finish? check its log")
        return
    off = max(plus1, minus1)
    msg = (f"{total:,} TC/AG conversion events (q>{quality}); at mask position: "
           f"{exact:,} ({exact / total:.1%}); at mask+1: {minus1:,}; at mask-1: {plus1:,}")
    if exact > 3 * max(off, 1):
        report("DYNAST", PASS, msg + " -> dynast genome_i and mask agree")
    elif off > 3 * max(exact, 1):
        report("DYNAST", FAIL, msg + " -> mask is shifted by 1 bp relative to dynast")
    else:
        report("DYNAST", WARN, msg + " -> few overlaps; check that this sample's donors are in the mask")


# ---------------------------------------------------------------------------
def tc_rates(run_dir):
    import pandas as pd
    path = os.path.join(run_dir, "rates.csv")
    if not os.path.exists(path):
        sys.exit(f"ERROR: {path} not found")
    df = pd.read_csv(path)
    return df.set_index("barcode")["TC"] if "barcode" in df else df["TC"]


def check_effect(unmasked, masked):
    a, b = tc_rates(unmasked), tc_rates(masked)
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]
    drop = 1 - b.median() / a.median() if a.median() > 0 else float("nan")
    lower = (b < a).mean()
    msg = (f"{len(common):,} cells; median T>C rate {a.median():.5f} -> {b.median():.5f} "
           f"({drop:.1%} lower); cells with lower rate after masking: {lower:.1%}")
    if (b > a + 1e-12).any():
        report("EFFECT", WARN, msg + "; some cells INCREASED (masking should never add conversions)")
    elif drop >= 0.01 and lower >= 0.9:
        report("EFFECT", PASS, msg)
    else:
        report("EFFECT", FAIL, msg + " -> masking barely changed anything; mask not applied?")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mask", required=True, help="SNP CSV passed to dynast --snp-csv")
    ap.add_argument("--fasta", help="genome FASTA (with .fai) used for alignment")
    ap.add_argument("--bam", help="BAM given to dynast (contig-name check)")
    ap.add_argument("--conversions", help="conversions.csv from an UNMASKED dynast count run")
    ap.add_argument("--unmasked", help="dynast count output dir WITHOUT --snp-csv")
    ap.add_argument("--masked", help="dynast count output dir WITH --snp-csv")
    ap.add_argument("--n-sample", type=int, default=20000, help="mask rows to test against FASTA")
    ap.add_argument("--quality", type=int, default=27, help="dynast count --quality (default 27)")
    ap.add_argument("--max-lines", type=int, default=0,
                    help="read only this many lines of conversions.csv (0 = all)")
    a = ap.parse_args()

    header, mask, n, problems = load_mask(a.mask)
    check_format(header, mask, n, problems, a.bam)
    if a.fasta:
        check_coordinates(mask, a.fasta, a.n_sample)
    if a.conversions:
        check_dynast(mask, a.conversions, a.quality, a.max_lines)
    if a.unmasked and a.masked:
        check_effect(a.unmasked, a.masked)

    n_fail = sum(s == FAIL for _, s in results)
    print(f"\nSummary: {len(results) - n_fail}/{len(results)} checks without FAIL")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
