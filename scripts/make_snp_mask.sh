#!/usr/bin/env bash
# =============================================================================
# make_snp_mask.sh
# Build dynast-ready SNP masks (genomic T>C and A>G) from per-donor WGS VCFs.
#
# Output format (what dynast's read_snp_csv expects for a 3-column file):
#     contig,genome_i,conversion        <- header line is REQUIRED
#     chr1,14652,TC                     <- genome_i is 0-BASED (= VCF POS - 1)
#
# Why 0-based: dynast records genome_i from pysam get_aligned_pairs(), which is
# 0-based, and masks a conversion only if genome_i matches exactly.
#
# Usage:  nohup bash make_snp_mask.sh > make_snp_mask.log 2>&1 &
# =============================================================================
set -euo pipefail

vcf_dir="${vcf_dir:-/path/to/wgs_vcfs}"
ref_fa="${ref_fa:-/path/to/reference/GRCh38.genome.fa}"
out_dir="${out_dir:-/path/to/WGS_SNP_mask_0-based}"

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

check_tools bcftools
chroms=$(echo chr{1..22} chrX chrY chrM | tr ' ' ',')
[[ -e "${ref_fa}.fai" ]] || { echo "ERROR: ${ref_fa}.fai missing (run samtools faidx)"; exit 1; }

vcf_dir=$(cd "${vcf_dir}" && pwd)              # absolute paths before cd
ref_fa=$(cd "$(dirname "${ref_fa}")" && pwd)/$(basename "${ref_fa}")
mkdir -p "${out_dir}"
cd "${out_dir}"

shopt -s nullglob
vcfs=("${vcf_dir}"/*.vcf.gz)
(( ${#vcfs[@]} )) || { echo "ERROR: no *.vcf.gz in ${vcf_dir}"; exit 1; }

for donor_vcf in "${vcfs[@]}"; do
    donor_id=$(basename "${donor_vcf}" .hard-filtered.vcf.gz)
    donor_out="${donor_id}_snps.csv"
    ropt=$(vcf_region_opt "${donor_vcf}")
    echo "Processing ${donor_id} ... (${ropt} region access)"

    # 1. primary chromosomes only
    # 2. split multiallelic records, left-normalize against the reference
    # 3. keep SNPs where THIS donor carries the ALT allele (0/1, 1/1);
    #    no FILTER restriction (for a mask, over-masking is cheap)
    # 4. keep genomic T>C and A>G, convert VCF POS (1-based) -> 0-based
    { echo "contig,genome_i,conversion"
      bcftools view "${ropt}" "${chroms}" "${donor_vcf}" -Ou \
      | bcftools norm -m -both -f "${ref_fa}" --check-ref s -Ou 2>/dev/null \
      | bcftools view -v snps -i 'GT="alt"' -Ou \
      | bcftools query -f '%CHROM,%POS,%REF,%ALT\n' \
      | awk -F',' 'BEGIN{OFS=","} ($3$4=="TC" || $3$4=="AG"){print $1,$2-1,$3$4}'
    } > "${donor_out}.tmp"
    mv "${donor_out}.tmp" "${donor_out}"

    n=$(( $(wc -l < "${donor_out}") - 1 ))
    echo "  ${n} T>C/A>G sites  (expect ~1.3-1.7 million per donor)"
done

# Union of all donors (skip each file's header before sorting)
union_out="snps_union_all.csv"
echo "contig,genome_i,conversion" > "${union_out}.tmp"
for f in ./*_snps.csv; do tail -n +2 "$f"; done | LC_ALL=C sort -u >> "${union_out}.tmp"
mv "${union_out}.tmp" "${union_out}"
echo "Union: $(( $(wc -l < "${union_out}") - 1 )) sites -> ${out_dir}/${union_out}"
