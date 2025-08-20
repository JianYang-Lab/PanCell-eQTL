#!/bin/bash
#SBATCH --qos=normal
#SBATCH --mem=20G
#SBATCH -c 4
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --output=LOG/merge_pca_%j.out
#SBATCH --error=LOG/merge_pca_%j.err
#SBATCH -J merge_pca

# ==============================================================================
# BED File Merging, Reference Intersection, and PCA Script
# ==============================================================================
# This script merges per-chromosome filtered PLINK BED files, intersects the
# resulting SNPs with a reference panel, and performs GRM/PCA analysis using GCTA.
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
module load plink/1.90 # Specific version for merging
module load plink/2.0  # For other tasks
module load bcftools

# TODO: Update this path to point to your GCTA executable.
GCTA_EXEC="/path/to/your/gcta-1.94.1"
NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Define base paths using the tissue argument
# TODO: Update this base path to your project's root directory.
BASE_PROJECT_DIR="/path/to/${tissue}"
VCF_DIR="${BASE_PROJECT_DIR}/VCF_TopMedimputed/filtered"
BED_DIR="${BASE_PROJECT_DIR}/VCF_TopMedimputed/filtered" # Assuming output goes to the same dir

# TODO: Update this path to your reference panel directory.
REF_DIR="/path/to/reference_panel"

# Create output and log directories if they don't exist
mkdir -p "${BED_DIR}"
mkdir -p "LOG"

# --- Main Workflow ---

# Create a temporary local directory on the compute node for faster I/O
LOCAL_DIR="/data/$$_${tissue}_merge_pca"
echo "Creating temporary directory: ${LOCAL_DIR}"
mkdir -p "${LOCAL_DIR}"

# Define a cleanup function to be called on script exit
cleanup() {
    echo "--- Cleaning up temporary files ---"
    rm -rf "${LOCAL_DIR}"
    echo "Temporary directory ${LOCAL_DIR} removed."
}
trap cleanup EXIT

echo "--- Starting Merge and PCA for ${tissue} ---"

# 1. Stage Data
echo "Step 1: Staging filtered BED files and reference files..."
cp "${VCF_DIR}"/merged_chr*_rsq03_mind01_msite01_maf01_hwe1e6.* "${LOCAL_DIR}/"
cp "${REF_DIR}"/ref.* "${LOCAL_DIR}/"
cd "${LOCAL_DIR}"

# 2. Merge per-chromosome BED files
echo "Step 2: Merging per-chromosome BED files..."
# Create a list of files to merge (chr2 through chr22)
for i in {2..22}; do
    echo "merged_chr${i}_rsq03_mind01_msite01_maf01_hwe1e6" >> mergelist.txt
done

# Use plink 1.9 for merging
plink --bfile merged_chr1_rsq03_mind01_msite01_maf01_hwe1e6 \
    --merge-list mergelist.txt \
    --make-bed --out merged_rsq03_mind01_msite01_maf01_hwe1e6 \
    --threads "${NUM_THREADS}"

# 3. Intersect with Reference Panel
echo "Step 3: Intersecting SNPs with reference panel..."
module load plink/2.0 # Switch back to plink2
plink2 --bfile merged_rsq03_mind01_msite01_maf01_hwe1e6 --write-snplist --out merged_rsq03_mind01_msite01_maf01_hwe1e6 --threads "${NUM_THREADS}"

plink2 --bfile merged_rsq03_mind01_msite01_maf01_hwe1e6 \
    --extract-intersect merged_rsq03_mind01_msite01_maf01_hwe1e6.snplist ref.snplist \
    --make-bed --out merged_isec --threads "${NUM_THREADS}"

plink2 --bfile ref \
    --extract-intersect merged_rsq03_mind01_msite01_maf01_hwe1e6.snplist ref.snplist \
    --make-bed --out ref_isec --threads "${NUM_THREADS}"

# 4. Merge with Reference, Prune, and run PCA
echo "Step 4: Merging with reference for PCA..."
module load plink/1.90 # Switch back for final merge
plink --bfile merged_isec --bmerge ref_isec.bed ref_isec.bim ref_isec.fam \
    --make-bed --out merged_ref_isec --threads "${NUM_THREADS}"

echo "Step 5: Pruning merged dataset for LD..."
module load plink/2.0
plink2 --bfile merged_ref_isec --indep-pairwise 50 5 0.9 --out merged_ref_pruned --threads "${NUM_THREADS}"
plink2 --bfile merged_ref_isec --extract merged_ref_pruned.prune.in --make-bed --out merged_ref_pruned --threads "${NUM_THREADS}"

echo "Step 6: Calculating GRM and running PCA with GCTA..."
"${GCTA_EXEC}" --bfile merged_ref_pruned --make-grm --out merged_ref_pruned --threads "${NUM_THREADS}"
"${GCTA_EXEC}" --grm merged_ref_pruned --pca 20 --out merged_ref_pruned_pca --threads "${NUM_THREADS}"

# 5. Copy Final Results Back
echo "Step 7: Copying final BED and PCA results to ${BED_DIR}..."
cp -v merged_ref_isec* "${BED_DIR}/"
cp -v merged_ref_pruned* "${BED_DIR}/"
cp -v merged_rsq03_mind01_msite01_maf01_hwe1e6.bed "${BED_DIR}/"
cp -v merged_rsq03_mind01_msite01_maf01_hwe1e6.bim "${BED_DIR}/"
cp -v merged_rsq03_mind01_msite01_maf01_hwe1e6.fam "${BED_DIR}/"

echo "--- Session completed successfully! ---"
