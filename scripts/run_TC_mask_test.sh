#!/usr/bin/env bash
# =============================================================================
# run_TC_mask_test.sh
# One-donor test of WGS SNP masking for dynast T>C counting:
#   old mask (raw VCF POS, 1-based)  vs  new mask (0-based)  vs  no mask
#
# All inputs are COPIED and all outputs are written to:
#   /path/to/TC_mask_test
#
# Usage:
#   nohup bash run_TC_mask_test.sh > run_TC_mask_test.log 2>&1 &
#   (activates the conda env "dynast" itself; needs dynast, bcftools,
#    samtools, pysam, pandas, numpy, matplotlib in that env)
# Re-running skips steps whose output already exists.
# =============================================================================
set -euo pipefail

# ------------------------------- settings ------------------------------------
donor="${donor:-DONOR}"
threads="${threads:-8}"

sample="${sample:-SAMPLE_day0}"
src_consensus_dir="${src_consensus_dir:-/path/to/${sample}/split_bam_consensus/${donor}}"
# An earlier masked dynast count in the consensus folder (if present). Only its
# small files (rates.csv, json/logs) are copied; it is added to the comparison plot.
src_existing_count="${src_existing_count:-${src_consensus_dir}/counts_SNP_corrected}"
src_vcf="${src_vcf:-/path/to/wgs_vcfs/${donor}.hard-filtered.vcf.gz}"
src_old_mask="${src_old_mask:-/path/to/WGS_SNP_mask_old_1based/${donor}_snps.csv}"
# Corrected mask already made by make_snp_mask.sh. If this file does not exist,
# the script builds the 0-based mask from the VCF instead (step 1b).
src_new_mask="${src_new_mask:-/path/to/WGS_SNP_mask_0-based/${donor}_snps.csv}"
src_fa="${src_fa:-/path/to/reference/GRCh38.genome.fa}"
src_gtf="${src_gtf:-/path/to/reference/gencode.annotation.gtf}"

# Pre-consensus donor BAM (split by donor, before dynast consensus), only used
# for the per-UMI plot. Leave empty to skip that plot.
src_align_bam="${src_align_bam:-}"

# dynast --tmp (large intermediate files while parsing the BAM). Must be on a
# disk with plenty of space -- NOT /tmp. Python's own socket files still go to
# TMPDIR=/tmp, which keeps those paths short (avoids "AF_UNIX path too long").
dynast_tmp="${dynast_tmp:-${base:-/path/to/TC_mask_test}/tmp}"
min_free_gb="${min_free_gb:-100}"     # warn if less than this is free

# Folder holding check_snp_mask.py, find_snp_examples.py, plot_mask_effect.py, plot_tc_reads.py
script_src="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

base="${base:-/path/to/TC_mask_test}"
# -----------------------------------------------------------------------------

in="${base}/00_input"; mk="${base}/01_mask"; cnt="${base}/02_dynast_count"
chk="${base}/03_checks"; plt="${base}/04_plots"; scr="${base}/scripts"; logs="${base}/logs"
mkdir -p "$in" "$mk" "$cnt" "$chk" "$plt" "$scr" "$logs"
export TMPDIR="/tmp"
step() { echo; echo "==== $(date '+%F %T')  $*"; }
cpi() { [[ -e "$2" ]] || cp "$1" "$2"; }   # copy if not there yet

# ------------------------- environment (conda: dynast) -----------------------
# Activates the conda env (default "dynast") if it is not already active, so the
# script also works under nohup / cron. Override: conda_env=myenv bash <script>
conda_env="${conda_env:-dynast}"
if [[ "${CONDA_DEFAULT_ENV:-}" != "$conda_env" ]]; then
    conda_base=$(conda info --base 2>/dev/null || true)
    for b in "$conda_base" "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" "$HOME/mambaforge"; do
        if [[ -n "$b" && -f "$b/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1091
            source "$b/etc/profile.d/conda.sh"; break
        fi
    done
    set +u; conda activate "$conda_env" || { echo "ERROR: cannot activate conda env '$conda_env'"; exit 1; }; set -u
fi

# bcftools / samtools: default = the ones in the env. If they fail (e.g.
# "libgsl.so.25: cannot open shared object file"), fix the env with
#   conda install -n dynast -c conda-forge -c bioconda bcftools samtools
# or point to another install: bcftools_bin=/path/bcftools samtools_bin=/path/samtools
bcftools_bin="${bcftools_bin:-bcftools}"
samtools_bin="${samtools_bin:-samtools}"
bcftools() { command "$bcftools_bin" "$@"; }
samtools() { command "$samtools_bin" "$@"; }

check_tools() {   # usage: check_tools bcftools samtools dynast python:pysam,pandas
    local t
    for t in "$@"; do
        if [[ "$t" == python:* ]]; then
            python -c "import ${t#python:}" 2>/dev/null \
                || { echo "ERROR: python cannot import ${t#python:} in env '$conda_env'"; exit 1; }
        elif [[ "$t" == dynast ]]; then
            # dynast has no --version flag (and "dynast -h" exits 1), so check
            # that the command exists and its Python package imports
            command -v dynast > /dev/null && python -c "import dynast" 2>/dev/null || {
                echo "ERROR: 'dynast' not found or not importable in env '$conda_env':"
                python -c "import dynast" 2>&1 | tail -2 || true; exit 1; }
        elif ! "$t" --version > /dev/null 2>&1; then
            echo "ERROR: '$t' does not run in env '$conda_env':"
            "$t" --version 2>&1 | head -3 || true
            exit 1
        fi
    done
    echo "env: ${CONDA_DEFAULT_ENV:-none}"
    for t in "$@"; do
        case "$t" in
            bcftools|samtools) echo "  $t: $("$t" --version | head -1)";;
            dynast)            echo "  dynast: $(python -c 'import dynast; print(dynast.__version__)' 2>/dev/null || echo '?')  ($(command -v dynast))";;
        esac
    done
}

# Region access needs a VCF index that is newer than the VCF. If the index is
# missing or older (htslib: "index file is older than the data file"), stream
# the file with -t instead of -r (slower, same result).
vcf_region_opt() {   # usage: vcf_region_opt file.vcf.gz  -> prints "-r" or "-t"
    local v="$1" i
    for i in "$v.tbi" "$v.csi"; do
        if [[ -e "$i" && ! "$i" -ot "$v" ]]; then echo "-r"; return; fi
    done
    echo "-t"
}
# -----------------------------------------------------------------------------

check_tools bcftools samtools dynast python:pysam,pandas,numpy,matplotlib

# ---------------------------------------------------------------- 0. inputs
step "0. Copy inputs to ${in}"
for f in "${src_consensus_dir}/consensus.bam" "$src_vcf" "$src_old_mask" "$src_fa" "$src_gtf"; do
    [[ -e "$f" ]] || { echo "ERROR: missing input $f"; exit 1; }
done
cpi "${src_consensus_dir}/consensus.bam" "${in}/consensus.bam"
if [[ -e "${src_consensus_dir}/consensus.bam.bai" ]]; then cpi "${src_consensus_dir}/consensus.bam.bai" "${in}/consensus.bam.bai"; fi
cp "${src_consensus_dir}"/run_info_*.json "${in}/" 2>/dev/null || true
cpi "$src_vcf" "${in}/$(basename "$src_vcf")"
cpi "$src_old_mask" "${in}/${donor}_snps_OLD_1based.csv"
if [[ -s "$src_new_mask" ]]; then cpi "$src_new_mask" "${in}/${donor}_snps_NEW_0based.csv"; fi
cpi "$src_fa" "${in}/genome.fa"
cpi "$src_gtf" "${in}/annotation.gtf"
if [[ -n "$src_align_bam" ]]; then cpi "$src_align_bam" "${in}/${donor}_align.bam"; fi
if [[ -s "${src_existing_count}/rates.csv" ]]; then
    ex="${in}/existing_counts_SNP_corrected"; mkdir -p "$ex"
    cp "${src_existing_count}/rates.csv" "$ex/"
    cp "${src_existing_count}"/*.json "${src_existing_count}"/*.log "$ex/" 2>/dev/null || true
    echo "  copied existing ${src_existing_count} (rates.csv + json/logs)"
fi
for s in check_snp_mask.py find_snp_examples.py plot_mask_effect.py plot_tc_reads.py run_TC_mask_test.sh; do
    [[ -e "${script_src}/${s}" ]] || { echo "ERROR: ${s} not found next to this script (${script_src})"; exit 1; }
    # skip if the script is already being run from the scripts/ folder
    [[ "${script_src}/${s}" -ef "${scr}/${s}" ]] || cp -f "${script_src}/${s}" "${scr}/"
done

bam="${in}/consensus.bam"
vcf="${in}/$(basename "$src_vcf")"
fa="${in}/genome.fa"
gtf="${in}/annotation.gtf"
old_raw="${in}/${donor}_snps_OLD_1based.csv"

# indexes (fresh copies -> also fixes "index file is older than the data file")
[[ -e "${bam}.bai" ]] || samtools index "$bam"
bcftools index -t -f "$vcf"
[[ -e "${fa}.fai" ]] || samtools faidx "$fa"
if [[ -n "$src_align_bam" && ! -e "${in}/${donor}_align.bam.bai" ]]; then samtools index "${in}/${donor}_align.bam"; fi

# show the original commands (check that GTF/strand match; see which --snp-csv
# the existing counts_SNP_corrected run used)
for j in "${in}"/run_info_*.json "${in}"/existing_counts_SNP_corrected/*.json; do
    [[ -e "$j" ]] && python -c "import json,sys; d=json.load(open(sys.argv[1])); c=d.get('call') or d.get('command') or d; print(sys.argv[1].split('/')[-2] + ':', str(c)[:600])" "$j" || true
done

# ---------------------------------------------------------------- 1. masks
step "1a. Old mask: add header only (positions stay 1-based, as dynast would read them)"
old_mask="${mk}/${donor}_snps_OLD_1based_with_header.csv"
if [[ "$(head -1 "$old_raw")" == "contig,genome_i,conversion" ]]; then
    cp "$old_raw" "$old_mask"
else
    { echo "contig,genome_i,conversion"; cat "$old_raw"; } > "$old_mask"
fi
echo "  $(( $(wc -l < "$old_mask") - 1 )) sites"

step "1b. New mask (0-based)"
new_mask="${mk}/${donor}_snps_0based.csv"
if [[ -s "${in}/${donor}_snps_NEW_0based.csv" ]]; then
    echo "  using existing mask: ${src_new_mask}"
    cp "${in}/${donor}_snps_NEW_0based.csv" "$new_mask"
elif [[ ! -s "$new_mask" ]]; then
    echo "  ${src_new_mask} not found -> building it from the VCF"
    chroms=$(echo chr{1..22} chrX chrY chrM | tr ' ' ',')
    { echo "contig,genome_i,conversion"
      bcftools view "$(vcf_region_opt "$vcf")" "$chroms" "$vcf" -Ou \
      | bcftools norm -m -both -f "$fa" --check-ref s -Ou \
      | bcftools view -v snps -i 'GT="alt"' -Ou \
      | bcftools query -f '%CHROM,%POS,%REF,%ALT\n' \
      | awk -F',' 'BEGIN{OFS=","} ($3$4=="TC" || $3$4=="AG"){print $1,$2-1,$3$4}'
    } > "$new_mask"
fi
echo "  $(( $(wc -l < "$new_mask") - 1 )) sites  (expect ~1.3-1.7 million)"

# ---------------------------------------------------------------- 2. dynast count
run_count() {   # $1 = output name, rest = extra dynast args
    local name="$1"; shift
    local out="${cnt}/${name}"
    if [[ -s "${out}/rates.csv" ]] && ! grep -q "An exception occurred" "${logs}/count_${name}.log" 2>/dev/null; then
        echo "  ${name}: done already, skipping"; return
    fi
    rm -rf "$out"; mkdir -p "$out"           # remove partial output of a failed run
    # The BAM parsing step is identical for all runs. dynast skips it when these
    # 4 files exist, so hard-link them from the unmasked run (no extra disk/time).
    local f parsed=(conversions.csv conversions.idx alignments.csv genes.pkl.gz)
    if [[ "$name" != TC ]]; then
        for f in "${parsed[@]}"; do ln -f "${cnt}/TC/${f}" "${out}/${f}" 2>/dev/null || cp "${cnt}/TC/${f}" "${out}/${f}"; done
    fi
    local tmpd="${dynast_tmp}/${name}"
    rm -rf "$tmpd"; mkdir -p "$dynast_tmp"   # dynast refuses to start if --tmp exists
    echo "  ${name}: running (log: ${logs}/count_${name}.log)"
    dynast count -t "$threads" --strand forward -g "$gtf" \
        --barcode-tag CB --umi-tag UB --conversion TC "$@" \
        -o "$out" --tmp "$tmpd" \
        "$bam" &> "${logs}/count_${name}.log" || true
    rm -rf "$tmpd"
    # dynast catches its own errors and still exits 0, so check the result
    if [[ ! -s "${out}/rates.csv" ]] || grep -q "An exception occurred" "${logs}/count_${name}.log"; then
        echo "ERROR: dynast count '${name}' failed. Last lines of ${logs}/count_${name}.log:"
        tail -n 40 "${logs}/count_${name}.log"
        grep -q "No space left on device" "${logs}/count_${name}.log" && disk_report
        exit 1
    fi
    echo "  ${name}: finished"
}
disk_report() {
    echo "Disk space:"
    df -h "$base" "$dynast_tmp" /tmp 2>/dev/null | awk 'NR==1 || !seen[$0]++'
}
mkdir -p "$dynast_tmp"
free_gb=$(df -Pk "$base" | awk 'NR==2{print int($4/1024/1024)}')
tmp_gb=$(df -Pk "$dynast_tmp" | awk 'NR==2{print int($4/1024/1024)}')
bam_gb=$(du -k "$bam" | awk '{printf "%.1f", $1/1024/1024}')
echo "  free: ${free_gb} GB at ${base}, ${tmp_gb} GB at dynast --tmp ${dynast_tmp}; consensus.bam = ${bam_gb} GB"
if (( free_gb < min_free_gb || tmp_gb < min_free_gb )); then
    echo "  WARNING: less than ${min_free_gb} GB free; dynast count may fail with 'No space left on device'."
    echo "           Set dynast_tmp=/path/on/a/big/disk or free space first."
fi
step "2. dynast count on consensus.bam (mask only, no --snp-threshold, for a fair comparison)"
run_count TC
run_count TC_rm_bg_old1based --snp-csv "$old_mask"
run_count TC_rm_bg_0based    --snp-csv "$new_mask"

# ---------------------------------------------------------------- 3. checks
step "3a. Check OLD mask (expected: COORDINATE, DYNAST and EFFECT fail)"
python "${scr}/check_snp_mask.py" --mask "$old_mask" --fasta "$fa" --bam "$bam" \
    --conversions "${cnt}/TC/conversions.csv" \
    --unmasked "${cnt}/TC" --masked "${cnt}/TC_rm_bg_old1based" \
    | tee "${chk}/check_OLD_mask.txt" || true

step "3b. Check NEW mask (expected: all PASS)"
python "${scr}/check_snp_mask.py" --mask "$new_mask" --fasta "$fa" --bam "$bam" \
    --conversions "${cnt}/TC/conversions.csv" \
    --unmasked "${cnt}/TC" --masked "${cnt}/TC_rm_bg_0based" \
    | tee "${chk}/check_NEW_mask.txt" || true

# ---------------------------------------------------------------- 4. plots
step "4a. Per-cell T>C rate: no mask vs old vs new"
runs=( --run "No mask=${cnt}/TC"
       --run "Old mask (1-based)=${cnt}/TC_rm_bg_old1based"
       --run "New mask (0-based)=${cnt}/TC_rm_bg_0based" )
[[ -s "${in}/existing_counts_SNP_corrected/rates.csv" ]] && \
    runs+=( --run "Existing counts_SNP_corrected=${in}/existing_counts_SNP_corrected" )
python "${scr}/plot_mask_effect.py" "${runs[@]}" \
    --title "${donor} (${sample}, consensus.bam)" \
    --out "${plt}/mask_effect_${donor}.pdf" | tee "${chk}/mask_effect_${donor}.txt"

step "4b. Expressed SNP sites that were counted as T>C without a mask"
python "${scr}/find_snp_examples.py" --conversions "${cnt}/TC/conversions.csv" \
    --mask "$new_mask" --gtf "$gtf" --top 20 | tee "${chk}/snp_examples_${donor}.tsv"

if [[ -n "$src_align_bam" ]]; then
    step "4c. Per-UMI plots (old vs new mask) for the top 3 genes"
    awk -F'\t' 'NR>1 && $7!="" && $7!="intergenic" && $1 ~ /^[0-9]+$/ && !seen[$7]++ {print $7"\t"$8}' \
        "${chk}/snp_examples_${donor}.tsv" | head -3 |
    while IFS=$'\t' read -r gene win; do
        for which in OLD NEW; do
            m=$old_mask; [[ $which == NEW ]] && m=$new_mask
            python "${scr}/plot_tc_reads.py" --bam "${in}/${donor}_align.bam" --fasta "$fa" --gtf "$gtf" \
                --gene "$gene" --xlim $win --snp-csv "$m" --snp-coord 0based \
                --out "${plt}/${gene}_${which}_mask.pdf" \
                > "${chk}/plot_${gene}_${which}_mask.txt" 2>&1 \
                && echo "  ${plt}/${gene}_${which}_mask.pdf" \
                || echo "  ${gene} ${which}: failed, see ${chk}/plot_${gene}_${which}_mask.txt"
        done
    done
else
    echo; echo "(4c skipped: set src_align_bam to the pre-consensus donor BAM to make per-UMI plots)"
fi

step "Done. Results in ${base}"
echo "  checks : ${chk}/check_OLD_mask.txt  ${chk}/check_NEW_mask.txt"
echo "  plots  : ${plt}/"
