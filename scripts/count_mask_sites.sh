#!/usr/bin/env bash
# =============================================================================
# count_mask_sites.sh
# Count TC and AG sites in each donor SNP mask and in the union mask, and the
# percentage of the genome they cover.
#
# Usage:  bash count_mask_sites.sh [mask_dir] [genome.fa]
# Output: table on screen + mask_site_counts.tsv in mask_dir
#
# Percentages (primary chromosomes chr1-22, X, Y, M):
#   pct_genome : sites / all non-N bases
#   pct_T_A    : sites / T+A bases (a TC site is always a T, an AG site an A,
#                so this is the fraction of possible T>C positions masked)
# The base counts are computed once from the FASTA (a few minutes) and cached
# in mask_dir/genome_base_counts.tsv.
# =============================================================================
set -euo pipefail

mask_dir="${1:-/path/to/WGS_SNP_mask_0-based}"
fasta="${2:-/path/to/reference/GRCh38.genome.fa}"
cd "$mask_dir"
out="mask_site_counts.tsv"
cache="genome_base_counts.tsv"

# ------------------------------------------------ genome base counts (cached)
if [[ ! -s "$cache" ]]; then
    echo "Counting bases in ${fasta} (primary chromosomes; once, cached in ${cache}) ..."
    python - "$fasta" > "${cache}.tmp" <<'EOF'
import sys
import pysam
fa = pysam.FastaFile(sys.argv[1])
chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
tot = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "length": 0}
for c in chroms:
    if c not in fa.references:
        continue
    s = fa.fetch(c).upper().encode()
    tot["length"] += len(s)
    for b in "ACGTN":
        tot[b] += s.count(b.encode())
non_n = tot["length"] - tot["N"]
print("length\tnon_N\tT_plus_A")
print(f"{tot['length']}\t{non_n}\t{tot['T'] + tot['A']}")
EOF
    mv "${cache}.tmp" "$cache"
fi
read -r g_len g_nonN g_TA < <(tail -n 1 "$cache")

# ------------------------------------------------ count sites
count_file() {   # prints: TC AG total  (skips the header line)
    awk -F',' 'NR>1 { n[$3]++; t++ } END { printf "%d\t%d\t%d\n", n["TC"]+0, n["AG"]+0, t+0 }' "$1"
}
row() {          # row <name> <TC> <AG> <total>
    awk -v n="$1" -v tc="$2" -v ag="$3" -v t="$4" -v g="$g_nonN" -v ta="$g_TA" \
        'BEGIN { printf "%s\t%d\t%d\t%d\t%.3f%%\t%.3f%%\n", n, tc, ag, t, 100*t/g, 100*t/ta }'
}

{
    printf "mask\tTC\tAG\ttotal\tpct_genome\tpct_T_A\n"
    sum_tc=0; sum_ag=0; sum_all=0; n_donor=0
    for f in *_snps.csv; do
        [[ "$f" == snps_union_all.csv ]] && continue
        read -r tc ag tot < <(count_file "$f")
        row "${f%_snps.csv}" "$tc" "$ag" "$tot"
        sum_tc=$((sum_tc + tc)); sum_ag=$((sum_ag + ag)); sum_all=$((sum_all + tot))
        n_donor=$((n_donor + 1))
    done
    if (( n_donor > 0 )); then
        row "MEAN_per_donor(n=${n_donor})" $((sum_tc / n_donor)) $((sum_ag / n_donor)) $((sum_all / n_donor))
        row "SUM_all_donors" "$sum_tc" "$sum_ag" "$sum_all"
    fi
    if [[ -e snps_union_all.csv ]]; then
        read -r tc ag tot < <(count_file snps_union_all.csv)
        row "UNION" "$tc" "$ag" "$tot"
    fi
} | tee "$out" | { if command -v column > /dev/null; then column -t; else cat; fi; }

echo
echo "Genome (primary chromosomes): ${g_len} bp, ${g_nonN} non-N bases, ${g_TA} T+A bases"
echo "Saved: ${mask_dir}/${out}"
