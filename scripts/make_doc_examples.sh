#!/usr/bin/env bash
# =============================================================================
# make_doc_examples.sh -- example figures for the Word guide, plot types 1 and 2.
# Run from the TC_mask_test folder:
#   bash scripts/make_doc_examples.sh        # plot types 1 and 2
#   bash scripts/make_doc_examples.sh 1      # only plot type 1 (basic)
#   bash scripts/make_doc_examples.sh 2      # only plot type 2 (mask before/after)
# Output: 04_plots/doc_examples/*.pdf  (+ one .txt log per figure)
# =============================================================================
set -uo pipefail
bam="${bam:-/path/to/split_bams/DONOR.bam}"
cons="${cons:-00_input/consensus.bam}"
mask="${mask:-01_mask/DONOR_snps_0based.csv}"
ad_m="${ad_m:-02_dynast_count/TC_rm_bg_0based/adata.h5ad}"
ad_u="${ad_u:-02_dynast_count/TC/adata.h5ad}"
out="${out:-04_plots/doc_examples}"
py="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/plot_tc_panel.py"
common=(--fasta 00_input/genome.fa --gtf 00_input/annotation.gtf --snp-csv "$mask"
        --format pdf --outdir "$out")
which="${1:-all}"
mkdir -p "$out"

run() {   # run <log name> <options...>
    local name="$1"; shift
    echo "== $name"
    python "$py" "$@" > "$out/$name.txt" 2>&1
    grep "^Saved:" "$out/$name.txt" || echo "   failed, see $out/$name.txt"
}

# ---- plot type 1: basic
if [[ "$which" == all || "$which" == 1 ]]; then
run 1a_basic_manual        --bam "$bam" "${common[@]}" --genes SOX2 GAPDH ACTB \
                           --vmax 4 --max-umis 40 --order position
run 1b_basic_allreads      --bam "$bam" "${common[@]}" --genes SOX2 GAPDH ACTB \
                           --vmax 4 --max-umis 40 --order position --reads-in-rows --row-height 0.15
run 1c_basic_ntr_auto      --bam "$bam" "${common[@]}" --adata "$ad_m" --no-mt-ribo \
                           --high-ntr-genes 1 --low-ntr-genes 1 --vmax 4 --max-umis 40
run 1d_basic_top_cells     --bam "$bam" "${common[@]}" --genes SOX2 GAPDH ACTB --top-cells 3 \
                           --vmax 4 --max-umis 40 --order position

fi

# ---- plot type 2: mask before/after (consensus.bam = what dynast counted)
if [[ "$which" == all || "$which" == 2 ]]; then
run 2a_mask_auto_genes     --bam "$cons" --min-bq 28 "${common[@]}" --adata "$ad_m" \
                           --adata-unmasked "$ad_u" --mask-genes 2 --no-mt-ribo \
                           --control-shift 1 --window 50 --vmax 4
run 2b_mask_GSTP1_top      --bam "$cons" --min-bq 28 "${common[@]}" --genes GSTP1 --before-after \
                           --top-cells 2 --control-shift 1 --window 50 --vmax 4
run 2c_mask_GSTP1_allreads --bam "$bam" "${common[@]}" --genes GSTP1 --before-after \
                           --control-shift 1 --window 50 --vmax 4 --reads-in-rows --row-height 0.25
run 2d_mask_only           --bam "$cons" --min-bq 28 "${common[@]}" --genes GSTP1 --before-after \
                           --only-masked --window 50 --vmax 4

fi

echo; ls -1 "$out"/*.pdf 2>/dev/null
