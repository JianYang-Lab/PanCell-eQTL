#!/bin/bash
#SBATCH -a 1-22
#SBATCH --qos=huge
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --output=LOG/vcfFilter_chr%a_%j.out
#SBATCH --error=LOG/vcfFilter_chr%a_%j.err
#SBATCH -J vcfFilter
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short

# ==============================================================================
# VCF Filtering Script
# ==============================================================================
# This script filters a merged, imputed VCF file for a specific chromosome.
# It filters based on imputation R-squared, MAF, HWE, and missingness.
# Designed to be run as a Slurm job array, one job per chromosome.
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error.
set -o pipefail # The exit status of a pipeline is the last command to exit with a non-zero status.

# --- Function Definitions ---

# Prints usage information for the script.
usage() {
    echo "Usage: $0 --tissue <tissue>"
    echo ""
    echo "Options:"
    echo "  --tissue    Name of the tissue (e.g., Colon) to specify the data directory."
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
module load plink/2.0
module load bcftools

# Define the chromosome for this job array task
readonly CHR_ID=${SLURM_ARRAY_TASK_ID}
readonly NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Define base paths using the tissue argument
# TODO: Update this path to point to your project's VCF directory.
readonly VCF_DIR="/path/to/your/project/${tissue}/VCF_1KGimputed/VCF_merged"

# Check if input directory exists
if [ ! -d "${VCF_DIR}" ]; then
    echo "Error: VCF directory not found at ${VCF_DIR}" >&2
    exit 1
fi

# Create log directory if it doesn't exist
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

echo "--- Starting VCF filtering for ${tissue}, Chromosome ${CHR_ID} ---"

# 1. Stage Data
cd "${LOCAL_DIR}"
readonly INPUT_VCF="${VCF_DIR}/merged_chr${CHR_ID}.vcf.gz"
if [ ! -f "${INPUT_VCF}" ]; then
    echo "Error: Input VCF not found: ${INPUT_VCF}" >&2
    exit 1
fi
echo "Staging input VCF to local directory..."
cp "${INPUT_VCF}" .

# 2. Filter by Imputation R-squared (R2 > 0.3)
echo "Step 2: Filtering by imputation R-squared > 0.3..."
readonly RSQ_VCF="merged_chr${CHR_ID}_rsq03.vcf.gz"
bcftools view -i 'INFO/R2>0.3' "merged_chr${CHR_ID}.vcf.gz" \
    -Oz -o "${RSQ_VCF}" --threads "${NUM_THREADS}"

# 3. PLINK2 Filtering and BED file creation
echo "Step 3: Filtering with PLINK2 (MAF, HWE, missingness)..."
readonly PLINK_OUT_PREFIX="merged_chr${CHR_ID}_filtered_rsq03"
plink2 --vcf "${RSQ_VCF}" \
    --set-all-var-ids '@:#:$r:$a' --rm-dup exclude-all 'list' --new-id-max-allele-len 1000 \
    --mind 0.01 --geno 0.01 --maf 0.01 --hwe 1e-6 \
    --double-id \
    --make-bed --threads "${NUM_THREADS}" \
    --out "${PLINK_OUT_PREFIX}"

# 4. Copy Results Back
echo "Step 4: Copying filtered PLINK files back to ${VCF_DIR}..."
cp -v "${PLINK_OUT_PREFIX}".* "${VCF_DIR}/"

echo "--- Session for chr${CHR_ID} completed successfully! ---"
