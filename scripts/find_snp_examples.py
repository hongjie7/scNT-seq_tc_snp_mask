#!/usr/bin/env python3
"""
find_snp_examples.py -- find expressed WGS SNP sites that dynast counts as T>C,
so you know which gene/position to plot with plot_tc_reads.py.

It tallies TC/AG conversion events from an UNMASKED dynast run (conversions.csv)
that fall exactly on a mask position (0-based, the corrected mask), and reports
the top sites with the gene they fall in. dynast stores plus-strand genomic bases, so a
T>C in the RNA is a genomic TC on a + strand gene or AG on a - strand gene. An
AG on a + strand gene (e.g. common mitochondrial A>G variants in MT-RNR1/2) is
an A>G in the RNA and does not affect the T>C rate; those are hidden unless --all.

Usage:
  python find_snp_examples.py \
      --conversions donor_count/<donor>/TC/conversions.csv \
      --mask WGS_SNP_mask_0-based/<donor>_snps.csv \
      --gtf gencode.v47.chr_patch_hapl_scaff.annotation.gtf --top 20
"""
import argparse
import bisect
import gzip
from collections import Counter, defaultdict


def load_mask(path):
    mask = defaultdict(lambda: defaultdict(set))
    with open(path) as f:
        for i, line in enumerate(f):
            p = line.rstrip("\n").split(",")
            if len(p) < 3 or not p[1].isdigit():
                continue
            mask[p[2]][p[0]].add(int(p[1]))
    return mask


def load_genes(gtf, chroms):
    genes = defaultdict(list)          # chrom -> [(start0, end, strand, name)]
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if p[2] != "gene" or p[0] not in chroms:
                continue
            attr = p[8]
            name = attr.split('gene_name "')[1].split('"')[0] if 'gene_name "' in attr else "?"
            gtype = attr.split('gene_type "')[1].split('"')[0] if 'gene_type "' in attr else ""
            genes[p[0]].append((int(p[3]) - 1, int(p[4]), p[6], name, gtype))
    for c in genes:
        genes[c].sort()
    return genes


def genes_at(genes, chrom, pos):
    """All genes overlapping pos: [(strand, name)], protein-coding first."""
    lst = genes.get(chrom, [])
    i = bisect.bisect_right(lst, (pos, float("inf")))
    hits = []
    for s, e, st, name, gtype in reversed(lst[max(0, i - 200):i]):
        if s <= pos < e:
            hits.append((gtype != "protein_coding", st, name))
    return [(st, n) for _, st, n in sorted(hits)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conversions", required=True, help="conversions.csv of an UNMASKED run")
    ap.add_argument("--mask", required=True, help="CORRECTED (0-based) mask CSV")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--quality", type=int, default=27)
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--max-lines", type=int, default=0, help="0 = read whole file")
    ap.add_argument("--all", action="store_true",
                    help="also list sites that do NOT count as T>C (A>G on a + strand gene, "
                         "T>C on a - strand gene, or intergenic)")
    a = ap.parse_args()

    mask = load_mask(a.mask)
    hits = Counter()
    with open(a.conversions) as f:
        for i, line in enumerate(f):
            if a.max_lines and i > a.max_lines:
                break
            p = line.rstrip("\n").split(",")
            if len(p) < 6 or p[4] not in ("TC", "AG") or not p[5].isdigit():
                continue
            if int(p[5]) <= a.quality:
                continue
            g = int(p[3])
            if g in mask.get(p[4], {}).get(p[2], ()):
                hits[(p[2], g, p[4])] += 1

    ranked = hits.most_common()
    genes = load_genes(a.gtf, {k[0] for k, _ in ranked[:5000]})
    print("rank\tchrom\tpos_1based\tgenome_i_0based\tconversion\tevents\tgene\tplot_window"
          "\tgene_strand\tcounted_as_TC")
    r = 0
    for (c, g, conv), n in ranked[:5000]:
        # dynast counts T>C as genomic TC on + strand genes and AG on - strand genes
        need = "+" if conv == "TC" else "-"
        ov = genes_at(genes, c, g)
        match = [nm for st, nm in ov if st == need]
        if match:
            gene, gstrand, counted = match[0], need, "yes"
        elif ov:
            gene, gstrand, counted = ov[0][1], ov[0][0], "no"
        else:
            gene, gstrand, counted = "intergenic", ".", "no"
        if counted == "no" and not a.all:
            continue
        r += 1
        print(f"{r}\t{c}\t{g + 1}\t{g}\t{conv}\t{n}\t{gene}\t{g + 1 - 50} {g + 1 + 50}"
              f"\t{gstrand}\t{counted}")
        if r >= a.top:
            break
    if not a.all:
        print("\n(only sites that dynast counts as T>C are listed; use --all to include "
              "A>G on + strand genes, T>C on - strand genes and intergenic sites)")
    print("\nPlot a site with:  python plot_tc_reads.py ... --gene <gene> --xlim <plot_window>")


if __name__ == "__main__":
    main()
