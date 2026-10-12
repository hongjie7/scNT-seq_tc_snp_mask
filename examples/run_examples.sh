#!/usr/bin/env bash
# Run from the tc_snp_mask folder:   bash examples/run_examples.sh
# Makes the simulated demo data, then one example of each plot type in examples/figures/.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
d="$here/demo_data"
py="$here/../scripts/plot_tc_panel.py"
[[ -e "$d/pre.bam" ]] || python "$here/make_demo_data.py"
common=(--fasta "$d/genome.fa" --gtf "$d/genes.gtf" --snp-csv "$d/mask_0based.csv"
        --vmax 4 --format png --outdir "$here/figures")

# 1 transcripts: a fast and a stable gene in the same cell
python "$py" --bam "$d/pre.bam" "${common[@]}" --genes FAST STABLE --max-umis 30
# 2 transcripts with every read in its row
python "$py" --bam "$d/pre.bam" "${common[@]}" --genes FAST --max-umis 15 --reads-in-rows --row-height 0.3
# 3 mask: no mask | corrected mask, zoomed on the SNP
python "$py" --bam "$d/consensus.bam" --min-bq 28 "${common[@]}" --genes FAST --before-after --window 60
# 4 scNT-seq layout
python "$py" --bam "$d/pre.bam" "${common[@]}" --genes FAST --scnt-layout --label-umis --max-umis 30
# 5 scNT-seq layout: no mask | mask
python "$py" --bam "$d/pre.bam" "${common[@]}" --genes FAST --scnt-layout --compare-layout mask \
    --window 60 --label-umis --max-umis 30
# 6 scNT-seq layout: pre-consensus | consensus.bam
python "$py" --bam "$d/pre.bam" "${common[@]}" --genes FAST --scnt-layout \
    --compare-consensus "$d/consensus.bam" --compare-layout side --label-umis --max-umis 30
# 7 transcripts: pre-consensus | consensus.bam
python "$py" --bam "$d/pre.bam" "${common[@]}" --genes FAST \
    --compare-consensus "$d/consensus.bam" --compare-panels masked --label-umis --max-umis 30

# short copies for the README
cd "$here/figures"
for pair in "basic_group-yes_mask-yes_FAST-STABLE:1_transcripts" "basic_allreads_group-yes_mask-yes_FAST:2_transcripts_reads" \
            "mask_before_after_group-yes_FAST:3_mask" "scnt_layout_group-yes_mask-yes_FAST:4_scnt_layout" \
            "inset_mask_group-yes_FAST:5_scnt_layout_mask" "scnt_layout_consensus_side_group-yes_mask-yes_FAST:6_scnt_layout_consensus" \
            "consensus_masked_group-yes_FAST:7_transcripts_consensus"; do
    f=$(ls "${pair%%:*}"_CELL*.png | head -1); mv "$f" "${pair##*:}.png"
done
ls -1
