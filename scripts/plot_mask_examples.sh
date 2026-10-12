#!/usr/bin/env bash
# =============================================================================
# plot_mask_examples.sh
# Before/after SNP-mask plots for a list of genes (one PDF per gene) plus a
# summary table of labeled molecules without mask, with the +1 bp shifted mask
# (negative control) and with the corrected 0-based mask.
#
# Usage (from the TC_mask_test folder):
#   bash scripts/plot_mask_examples.sh                      # default gene list
#   bash scripts/plot_mask_examples.sh GSTP1 PTMA TXN       # your genes
# Settings can be overridden: bam=... mask=... window=... bash scripts/plot_mask_examples.sh
#   panels=masked        only the corrected 0-based mask panel
#   panels=before_after  no mask + corrected mask (no negative control)
# =============================================================================
set -uo pipefail

bam="${bam:-00_input/consensus.bam}"                 # consensus = what dynast counted
mask="${mask:-01_mask/DONOR_snps_0based.csv}"
fasta="${fasta:-00_input/genome.fa}"
gtf="${gtf:-00_input/annotation.gtf}"
window="${window:-50}"
min_bq="${min_bq:-28}"                               # dynast counts quality > 27
max_umis="${max_umis:-60}"                           # transcripts per panel
subsample="${subsample:-labeled_first}"              # labeled_first | proportional | first
order="${order:-labeled_first}"                      # labeled_first | position
outdir="${outdir:-04_plots/mask_examples}"
# panels: "all" = no mask + shift control + corrected mask
#         "before_after" = no mask + corrected mask
#         "masked" = corrected (0-based) mask only
panels="${panels:-all}"
case "$panels" in
    all)          extra=(--control-shift 1) ;;
    before_after) extra=() ;;
    masked)       extra=(--only-masked) ;;
    *) echo "ERROR: panels must be all, before_after or masked"; exit 1 ;;
esac
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

genes=("$@")
(( ${#genes[@]} )) || genes=(GSTP1 PTMA PSMA7 PEBP1 RPS21 TXN)

mkdir -p "$outdir"
summary="${outdir}/summary.tsv"
printf "gene\tsite\tcell\tUMIs\tlabeled_no_mask\tlabeled_shift_control\tlabeled_with_mask\tTC_removed\tstatus\n" > "$summary"

for g in "${genes[@]}"; do
    log="${outdir}/${g}.txt"
    pdf="${outdir}/${panels}_${g}.pdf"
    python "${script_dir}/plot_tc_panel.py" --bam "$bam" --fasta "$fasta" --gtf "$gtf" \
        --genes "$g" --before-after --min-bq "$min_bq" --snp-csv "$mask" \
        "${extra[@]}" --window "$window" --max-umis "$max_umis" --subsample "$subsample" --order "$order" --out "$pdf" > "$log" 2>&1
    rc=$?
    saved=$(grep -m1 "^Saved:" "$log" | sed 's/^Saved: //')
    [[ -n "$saved" ]] && pdf="$saved"      # file name includes the cell barcode
    site=$(grep -m1 "^Showing" "$log" | awk '{print $2}')
    cell=$(grep -m1 "^Showing" "$log" | sed -E 's/.* in cell ([^,]+),.*/\1/')
    get() { grep -m1 "^$g $1" "$log" | sed -E 's/.*: ([0-9]+) UMIs, ([0-9]+) labeled.*, ([0-9]+) T>C events removed/\1 \2 \3/'; }
    read -r n l0 _ <<< "$(get 'no mask')"
    read -r _ ls _ <<< "$(get 'mask shifted')"
    read -r n1 l1 rm <<< "$(get 'new 0-based mask')"
    n="${n:-$n1}"
    if (( rc == 0 )) && grep -q "does not change" "$log"; then
        status="ok_mask_no_effect"
        echo "made ${pdf}   (${g}: the mask does not change this gene here; ${l1:-?} of ${n:-?} labeled)"
    elif (( rc == 0 )); then
        status="ok"
        msg="${l1:-?} of ${n:-?} labeled with mask, ${rm:-?} T>C removed"
        [[ -n "${l0:-}" ]] && msg="${l0} -> ${msg}"
        [[ -n "${ls:-}" ]] && msg="${msg}; shift control ${ls} labeled"
        echo "made ${pdf}   (${msg})"
    elif grep -q "does not change" "$log"; then
        status="no_masked_TC"
        echo "${g}: the mask does not change this gene in this sample (no plot)"
    else
        status="error"
        echo "${g}: failed, see ${log}"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$g" "${site:-}" "${cell:-}" "${n:-}" \
        "${l0:-}" "${ls:-}" "${l1:-}" "${rm:-}" "$status" >> "$summary"
done

echo
column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
echo "Summary: ${summary}"
