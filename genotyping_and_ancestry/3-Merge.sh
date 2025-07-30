#!/bin/bash
#SBATCH -a 1-22
#SBATCH --qos=huge
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --output=LOG/merge_chr%a_%j.out
#SBATCH --error=LOG/merge_chr%a_%j.err
#SBATCH -J merge_vcf
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short

# ==============================================================================
# VCF Merging Script
# ==============================================================================
# This script merges per-sample imputed VCF files for a specific chromosome
# using bcftools. It is designed to be run as a Slurm job array.
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error.
set -o pipefail # Exit status of a pipeline is the last command to exit with a non-zero status.

# --- Function Definitions ---

# Prints usage information for the script.
usage() {
    echo "Usage: $0 --tissue <tissue>"
    echo ""
    echo "Options:"
    echo "  --tissue    Name of the tissue (e.g., Colon) corresponding to the data directory."
    echo "  -h, --help  Display this help message."
    exit 1
}

# --- Configuration & Argument Parsing ---

# Parse command-line arguments
if [ "$#" -eq 0 ]; then
    usage
fi

while [ "$#" -gt 0 ]; do
    case "$1" in
        --tissue)   tissue="$2"; shift 2;;
        -h|--help)  usage;;
        *)          echo "Unknown option: $1" >&2; usage;;
    esac
done

# Check if the tissue argument was provided
if [ -z "${tissue:-}" ]; then
    echo "Error: --tissue argument is required." >&2
    usage
fi

# --- Environment and Path Setup ---

echo "--- Setting up environment ---"
module load bcftools

# Define the chromosome for this job array task
readonly CHR_ID=${SLURM_ARRAY_TASK_ID}
readonly NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Define base paths using the tissue argument
readonly BASE_DIR="/path/to/your/project/${tissue}"
readonly INPUT_DIR="${BASE_DIR}/VCF_1KGimputed/VCF_persample"
readonly OUTPUT_DIR="${BASE_DIR}/VCF_1KGimputed/VCF_merged"

# Check if input directory exists
if [ ! -d "${INPUT_DIR}" ]; then
    echo "Error: Input directory not found at ${INPUT_DIR}" >&2
    exit 1
fi

# Create output and log directories if they don't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "LOG"

# --- Main Workflow ---

# Create a temporary local directory on the compute node for faster I/O
readonly LOCAL_DIR="/data/$$_${tissue}_chr${CHR_ID}"
echo "Creating temporary directory: ${LOCAL_DIR}"
mkdir -p "${LOCAL_DIR}"

# Define a cleanup function to be called on script exit
cleanup() {
    echo "--- Cleaning up temporary files ---"
    rm -rf "${LOCAL_DIR}"
    echo "Temporary directory ${LOCAL_DIR} removed."
}
trap cleanup EXIT

echo "--- Starting VCF merge for ${tissue}, Chromosome ${CHR_ID} ---"

# 1. Generate File List and Stage Data
cd "${INPUT_DIR}"
readonly SAMPLE_LIST="sample_chr${CHR_ID}.list"

echo "Generating file list for chr${CHR_ID}..."
# Find all imputed VCF files for the current chromosome
ls *_chr${CHR_ID}_imputed.vcf.gz > "${LOCAL_DIR}/${SAMPLE_LIST}"

if [ ! -s "${LOCAL_DIR}/${SAMPLE_LIST}" ]; then
    echo "Warning: No imputed VCF files found for chr${CHR_ID}. Exiting."
    exit 0
fi

echo "Staging VCF files to ${LOCAL_DIR}..."
# Copy the relevant VCF and index files to the local directory
# Using rsync for potentially better performance with many small files
rsync -a --files-from="${LOCAL_DIR}/${SAMPLE_LIST}" . "${LOCAL_DIR}/"
# Also copy the corresponding index files
sed 's/\.gz$/.gz.csi/' "${LOCAL_DIR}/${SAMPLE_LIST}" | xargs -I {} rsync -a --ignore-missing-args {} "${LOCAL_DIR}/"

# 2. Merge VCFs
cd "${LOCAL_DIR}"
echo "Merging VCFs for chr${CHR_ID} using ${NUM_THREADS} threads..."
readonly MERGED_VCF="merged_chr${CHR_ID}.vcf.gz"

bcftools merge \
    -l "${SAMPLE_LIST}" \
    --threads "${NUM_THREADS}" \
    -Oz -o "${MERGED_VCF}"

echo "Indexing the merged VCF..."
bcftools index "${MERGED_VCF}"

echo "Merging finished."

# 3. Copy Results Back
echo "Copying merged VCF back to ${OUTPUT_DIR}..."
cp -v "${MERGED_VCF}"* "${OUTPUT_DIR}/"

echo "--- Session for chr${CHR_ID} completed successfully! ---"

# The 'trap cleanup EXIT' will handle the removal of the local directory.
