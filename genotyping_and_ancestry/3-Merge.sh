#!/bin/bash
#SBATCH -a 1-22
#SBATCH --qos=huge
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --output=LOG/merge_chr%a_%j.out
#SBATCH --error=LOG/merge_chr%a_%j.err
#SBATCH -J merge_vcf
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short

# Exit on error, undefined variable, or pipeline failure
set -euo pipefail

# --- 1. Argument and Environment Setup ---

# Simplified argument parsing for a single required flag
if [[ "${1:-}" != "--tissue" || -z "${2:-}" ]]; then
    echo "Error: Missing or invalid arguments." >&2
    echo "Usage: $0 --tissue <tissue_name>" >&2
    exit 1
fi
tissue="$2"

# Load necessary software
module load bcftools

# Define key variables from Slurm environment
CHR_ID=${SLURM_ARRAY_TASK_ID}
NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Define paths
INPUT_DIR="/path/to/${tissue}/VCF_1KGimputed/VCF_persample"
OUTPUT_DIR="/path/to/${tissue}/VCF_1KGimputed/VCF_merged"

# Ensure directories exist
mkdir -p "${OUTPUT_DIR}" LOG
if [ ! -d "${INPUT_DIR}" ]; then
    echo "Error: Input directory not found at ${INPUT_DIR}" >&2
    exit 1
fi

# --- 2. Temporary Workspace & Cleanup ---

# Use fast, node-local storage for temporary files to improve I/O speed.
LOCAL_DIR=$(mktemp -d "/data/merge_chr${CHR_ID}.XXXXXXXX")

# 'trap' ensures the cleanup function runs automatically when the script exits
cleanup() {
    echo "--- Cleaning up ${LOCAL_DIR} ---"
    rm -rf "${LOCAL_DIR}"
}
trap cleanup EXIT

# --- 3. Main Workflow ---

echo "--- Starting merge for ${tissue}, Chromosome ${CHR_ID} ---"

# Stage each VCF and its index. The list consumed by bcftools names the
# staged copies, so the merge reads only from node-local storage.
shopt -s nullglob
vcfs=("${INPUT_DIR}"/*_chr"${CHR_ID}"_imputed.vcf.gz)
if (( ${#vcfs[@]} == 0 )); then
    echo "Warning: No imputed VCFs found for chr${CHR_ID}. Exiting."
    exit 0
fi

SAMPLE_LIST="${LOCAL_DIR}/sample_chr${CHR_ID}.list"
for vcf in "${vcfs[@]}"; do
    staged_vcf="${LOCAL_DIR}/${vcf##*/}"
    cp -p -- "${vcf}" "${staged_vcf}"
    if [[ -s "${vcf}.csi" && ! "${vcf}" -nt "${vcf}.csi" ]]; then
        cp -p -- "${vcf}.csi" "${staged_vcf}.csi"
    elif [[ -s "${vcf}.tbi" && ! "${vcf}" -nt "${vcf}.tbi" ]]; then
        cp -p -- "${vcf}.tbi" "${staged_vcf}.tbi"
    else
        bcftools index --threads "${NUM_THREADS}" "${staged_vcf}"
    fi
    printf '%s\n' "${staged_vcf}" >> "${SAMPLE_LIST}"
done

MERGED_VCF="${LOCAL_DIR}/merged_chr${CHR_ID}.vcf.gz"
if (( ${#vcfs[@]} == 1 )); then
    cp -- "${staged_vcf}" "${MERGED_VCF}"
else
    echo "Merging VCFs using ${NUM_THREADS} threads..."
    bcftools merge \
        -l "${SAMPLE_LIST}" \
        --threads "${NUM_THREADS}" \
        -Oz -o "${MERGED_VCF}"
fi

echo "Indexing merged VCF..."
bcftools index --threads "${NUM_THREADS}" "${MERGED_VCF}"

# Copy final results back to the persistent storage directory
echo "Copying final results to ${OUTPUT_DIR}..."
cp "${MERGED_VCF}"* "${OUTPUT_DIR}/"

echo "--- Job for chr${CHR_ID} completed successfully! ---"