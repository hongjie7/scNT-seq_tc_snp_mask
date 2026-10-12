#!/usr/bin/env bash
# =============================================================================
# run_dynast_bg_correction.sh
# For each sample, run dynast count twice on the consensus BAM:
#   TC        : no SNP mask (reference for comparison)
#   TC_rm_bg  : WGS SNP mask (union of all donors, 0-based) + de novo SNP detection
# Then validate with check_snp_mask.py.
#
# Usage:  nohup bash run_dynast_bg_correction.sh > run_dynast_bg_correction.log 2>&1 &
# =============================================================================
set -euo pipefail

threads="${threads:-16}"
gtf="${gtf:-/path/to/reference/gencode.annotation.gtf}"
fa="${fa:-/path/to/reference/GRCh38.genome.fa}"
snp_csv="${snp_csv:-/path/to/WGS_SNP_mask_0-based/snps_union_all.csv}"
CONSENSUS_BASE="${CONSENSUS_BASE:-/path/to/dynast_results}"
OUTPUT_BASE="${OUTPUT_BASE:-/path/to/dynast_results/WGS_SNP_bg_correction}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TMPDIR="/tmp"   # only Python's small socket files go here (keeps paths short)
# dynast --tmp holds large intermediate files: keep it OFF /tmp (root disk is small)
dynast_tmp="${dynast_tmp:-${OUTPUT_BASE}/tmp}"

samples=(
  "SAMPLE_day0"
  "SAMPLE_day2"
  "SAMPLE_day4"
  "SAMPLE_day6"
  "SAMPLE_day8"
  "SAMPLE_day16"
)

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

check_tools dynast python:pysam,pandas
[[ -e "${script_dir}/check_snp_mask.py" ]] || { echo "ERROR: check_snp_mask.py not next to this script"; exit 1; }

# --- guard: dynast needs this exact header for a 3-column SNP CSV ----------
[[ "$(head -1 "${snp_csv}")" == "contig,genome_i,conversion" ]] || {
    echo "ERROR: ${snp_csv} must start with 'contig,genome_i,conversion'"; exit 1; }

run_count() {   # run_count <sample_out_dir> <name> [extra dynast args]
    local sdir="$1" name="$2"; shift 2
    local out="${sdir}/${name}" log="${sdir}/count_${name}.log"
    if [[ -s "${out}/rates.csv" ]] && ! grep -q "An exception occurred" "$log" 2>/dev/null; then
        echo "  ${name}: done already, skipping"; return
    fi
    rm -rf "$out"; mkdir -p "$out"
    local f
    if [[ "$name" != TC ]]; then
        for f in conversions.csv conversions.idx alignments.csv genes.pkl.gz; do
            ln -f "${sdir}/TC/${f}" "${out}/${f}" 2>/dev/null || cp "${sdir}/TC/${f}" "${out}/${f}"
        done
    fi
    local tmpd="${dynast_tmp}/${sample}_${name}"
    rm -rf "$tmpd"; mkdir -p "$dynast_tmp"
    echo "  ${name}: running (log: ${log})"
    dynast count "${common[@]}" "$@" -o "$out" --tmp "$tmpd" "${bam}" &> "$log" || true
    rm -rf "$tmpd"
    # dynast catches its own errors and exits 0, so check the result
    if [[ ! -s "${out}/rates.csv" ]] || grep -q "An exception occurred" "$log"; then
        echo "ERROR: dynast count ${sample}/${name} failed. Last lines of ${log}:"
        tail -n 30 "$log"
        df -h "$OUTPUT_BASE" "$dynast_tmp" 2>/dev/null
        exit 1
    fi
}

for sample in "${samples[@]}"; do
    echo ">>> $(date '+%F %T') ${sample}"
    bam="${CONSENSUS_BASE}/${sample}/consensus/consensus.bam"
    barcodes="${CONSENSUS_BASE}/${sample}/cells_with_greater_than_700_UMIs.txt"
    out="${OUTPUT_BASE}/${sample}"
    mkdir -p "${out}"

    common=( -t "${threads}" -g "${gtf}" --strand forward
             --barcode-tag CB --umi-tag UB --conversion TC
             --barcodes "${barcodes}" )

    # 1) unmasked
    run_count "${out}" TC

    # 2) WGS-masked (+ de novo detection of high-frequency sites).
    #    BAM parsing is identical, so reuse the parsed files of the unmasked run.
    run_count "${out}" TC_rm_bg \
        --snp-csv "${snp_csv}" --snp-threshold 0.5 --snp-min-coverage 5

    # 3) validate
    python "${script_dir}/check_snp_mask.py" \
        --mask "${snp_csv}" --fasta "${fa}" --bam "${bam}" \
        --conversions "${out}/TC/conversions.csv" \
        --unmasked "${out}/TC" --masked "${out}/TC_rm_bg" \
        | tee "${out}/check_snp_mask.txt" || true
done
echo "All samples processed."
