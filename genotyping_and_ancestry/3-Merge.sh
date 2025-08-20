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
if [[ "$1" != "--tissue" || -z "$2" ]]; then
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
# The directory is unique to this job ($$) and array task (${CHR_ID}).
LOCAL_DIR="/data/$$_${tissue}_chr${CHR_ID}"
mkdir -p "${LOCAL_DIR}"

# 'trap' ensures the cleanup function runs automatically when the script exits
cleanup() {
    echo "--- Cleaning up ${LOCAL_DIR} ---"
    rm -rf "${LOCAL_DIR}"
}
trap cleanup EXIT

# --- 3. Main Workflow ---

echo "--- Starting merge for ${tissue}, Chromosome ${CHR_ID} ---"

# Find all relevant VCF files and create a list inside the temporary directory
SAMPLE_LIST="${LOCAL_DIR}/sample_chr${CHR_ID}.list"
find "${INPUT_DIR}" -maxdepth 1 -name "*_chr${CHR_ID}_imputed.vcf.gz" > "${SAMPLE_LIST}"

# Exit gracefully if no VCF files were found for this chromosome
if [ ! -s "${SAMPLE_LIST}" ]; then
    echo "Warning: No imputed VCFs found for chr${CHR_ID}. Exiting."
    exit 0
fi

# Stage the VCFs and their corresponding index files to the local directory
rsync -a --files-from="${SAMPLE_LIST}" "${INPUT_DIR}/" "${LOCAL_DIR}/"
sed 's/\.gz$/.gz.csi/' "${SAMPLE_LIST}" | xargs -I {} rsync -a --ignore-missing-args "${INPUT_DIR}/"{} "${LOCAL_DIR}/"

# Perform the merge and index operations inside the fast local directory
MERGED_VCF="${LOCAL_DIR}/merged_chr${CHR_ID}.vcf.gz"

echo "Merging VCFs using ${NUM_THREADS} threads..."
bcftools merge \
    -l "${SAMPLE_LIST}" \
    --threads "${NUM_THREADS}" \
    -Oz -o "${MERGED_VCF}"

echo "Indexing merged VCF..."
bcftools index --threads "${NUM_THREADS}" "${MERGED_VCF}"

# Copy final results back to the persistent storage directory
echo "Copying final results to ${OUTPUT_DIR}..."
cp "${MERGED_VCF}"* "${OUTPUT_DIR}/"

echo "--- Job for chr${CHR_ID} completed successfully! ---"