#!/usr/bin/env bash
# =============================================================================
# plot_top_cells.sh
# For each gene, take the N cells with the most UMIs and make, per cell:
#   1) pre-consensus vs consensus.bam (two panels, same rows, highlighted UMI)
#   2) scNT-seq layout figure (--scnt-layout) for the same cell + UMI (inset: its reads + consensus record)
#
# Usage (from the TC_mask_test folder):
#   bash scripts/plot_top_cells.sh                      # SOX2 ACTB, top 5 cells
#   n=3 bash scripts/plot_top_cells.sh GAPDH JUN
# Settings: bam= cons= mask= n= panels= (masked|nomask) max_umis= row_height= outdir=
# =============================================================================
set -uo pipefail

bam="${bam:-/path/to/split_bams/DONOR.bam}"
cons="${cons:-00_input/consensus.bam}"
mask="${mask:-01_mask/DONOR_snps_0based.csv}"
fasta="${fasta:-00_input/genome.fa}"
gtf="${gtf:-00_input/annotation.gtf}"
n="${n:-5}"
panels="${panels:-masked}"
max_umis="${max_umis:-60}"
row_height="${row_height:-0.15}"
label="${label:-Day 0}"
outdir="${outdir:-04_plots/top_cells}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py="${script_dir}/plot_tc_panel.py"

genes=("$@")
(( ${#genes[@]} )) || genes=(SOX2 ACTB)
mkdir -p "$outdir"
summary="${outdir}/summary.tsv"
printf "gene\trank\tcell\tUMIs\thighlighted_UMI\tcompare_pdf\tscnt_layout_pdf\n" > "$summary"

for gene in "${genes[@]}"; do
    echo "===== $gene: top $n cells"
    mapfile -t tops < <(python "$py" --bam "$bam" --fasta "$fasta" --gtf "$gtf" \
                          --genes "$gene" --list-cells "$n" 2>/dev/null | grep "^TOPCELL")
    if (( ${#tops[@]} == 0 )); then
        echo "  no cells with $gene UMIs (gene name in the GTF?)"; continue
    fi
    rank=0
    for line in "${tops[@]}"; do
        rank=$((rank + 1))
        IFS=$'\t' read -r _ _ cell numi <<< "$line"
        base="${outdir}/${gene}_top${rank}"
        echo "--- $gene #$rank  cell $cell  ($numi UMIs)"

        python "$py" --bam "$bam" --compare-consensus "$cons" --compare-panels "$panels" \
            --fasta "$fasta" --gtf "$gtf" --genes "$gene" --cell "$cell" --snp-csv "$mask" \
            --label-umis --vmax 4 --max-umis "$max_umis" --row-height "$row_height" \
            --out "${base}_pre_vs_consensus.pdf" > "${base}_pre_vs_consensus.txt" 2>&1
        umi=$(grep -m1 "^Highlighted transcript: UMI " "${base}_pre_vs_consensus.txt" | awk '{print $4}')
        pdf1=$(grep -m1 "^Saved:" "${base}_pre_vs_consensus.txt" | sed 's/^Saved: //')

        python "$py" --bam "$bam" --compare-consensus "$cons" \
            --fasta "$fasta" --gtf "$gtf" --genes "$gene" --cell "$cell" \
            ${umi:+--umi-reads "$umi"} --scnt-layout --label-umis --sample-label "$label" \
            --snp-csv "$mask" --max-umis "$max_umis" --row-height "$row_height" \
            --out "${base}_scnt_layout.pdf" > "${base}_scnt_layout.txt" 2>&1
        pdf2=$(grep -m1 "^Saved:" "${base}_scnt_layout.txt" | sed 's/^Saved: //')

        echo "    UMI ${umi:-none}:  ${pdf1:-FAILED (see ${base}_pre_vs_consensus.txt)}"
        echo "                      ${pdf2:-FAILED (see ${base}_scnt_layout.txt)}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$gene" "$rank" "$cell" "$numi" "${umi:-}" \
            "${pdf1:-}" "${pdf2:-}" >> "$summary"
    done
done
echo
echo "Summary: $summary"
