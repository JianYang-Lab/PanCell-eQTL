#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=25G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH --qos=normal
#SBATCH --job-name=monopogen
#SBATCH --output=monopogen_%j.log

# ==============================================================================
# Genotyping using Monopogen
# ==============================================================================
# This script automates the process of running the Monopogen pipeline on
# single-cell BAM files to call germline SNVs. It handles data staging,
# preprocessing, and germline variant calling.
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error.
set -o pipefail # Exit status of a pipeline is the last command to exit with a non-zero status.

# --- Function Definitions ---

# Prints usage information for the script.
usage() {
    echo "Usage: $0 --tissue <tissue> --project <project> --sample <sample>"
    echo ""
    echo "Options:"
    echo "  --tissue    Name of the tissue (e.g., Colon)."
    echo "  --project   Name of the project (e.g., PRJNA123456)."
    echo "  --sample    Name of the sample (e.g., Sample1)."
    echo "  -h, --help  Display this help message."
    exit 1
}

# --- Configuration & Argument Parsing ---

# Static paths for software and reference data
readonly MONOPOGEN_PATH="/path/to/SOFTWARE/Monopogen"
readonly REFERENCE_PATH="/path/to/genome_reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
readonly IMPUTATION_PATH="/path/to/LD_reference/1000G/1KGP3_3202_GRCh38"

# Parse command-line arguments
if [ "$#" -eq 0 ]; then
    usage
fi

while [ "$#" -gt 0 ]; do
    case "$1" in
        --tissue)   tissue="$2"; shift 2;;
        --project)  project="$2"; shift 2;;
        --sample)   sample="$2"; shift 2;;
        -h|--help)  usage;;
        *)          echo "Unknown option: $1" >&2; usage;;
    esac
done

# Check if all required arguments were provided
if [ -z "${tissue:-}" ] || [ -z "${project:-}" ] || [ -z "${sample:-}" ]; then
    echo "Error: Missing one or more required arguments." >&2
    usage
fi

# Construct project path and check for existence
readonly PROJECT_PATH="/path/to/your/project/${tissue}/raw_data/${project}"
if [ ! -d "${PROJECT_PATH}" ]; then
    echo "Error: Project path does not exist: ${PROJECT_PATH}" >&2
    exit 1
fi

# --- Environment Setup ---

echo "--- Setting up environment ---"
module load samtools
export LD_LIBRARY_PATH="${MONOPOGEN_PATH}/apps:${LD_LIBRARY_PATH}"
echo "LD_LIBRARY_PATH set."

# --- Main Workflow ---

# Create a temporary local directory on the compute node for faster I/O
readonly LOCAL_DIR="/data/$$_${sample}"
echo "Creating temporary directory: ${LOCAL_DIR}"
mkdir -p "${LOCAL_DIR}"

# Define a cleanup function to be called on script exit
cleanup() {
    echo "--- Cleaning up temporary files ---"
    rm -rf "${LOCAL_DIR}"
    echo "Temporary directory ${LOCAL_DIR} removed."
}
trap cleanup EXIT

# 1. Stage Data to Local Directory
echo "--- 1. Staging BAM file to local directory ---"
BAM_NAME="${project}_${sample}.bam"
BAM_FOUND=false

# Check for BAM file in a few common locations
if [ -f "${PROJECT_PATH}/bam/${BAM_NAME}" ]; then
    cp "${PROJECT_PATH}/bam/${BAM_NAME}"* "${LOCAL_DIR}/"
    BAM_FOUND=true
elif [ -f "${PROJECT_PATH}/${sample}/outs/${BAM_NAME}" ]; then
    cp "${PROJECT_PATH}/${sample}/outs/${BAM_NAME}"* "${LOCAL_DIR}/"
    BAM_FOUND=true
elif [ -f "${PROJECT_PATH}/${sample}/outs/possorted_genome_bam.bam" ]; then
    cp "${PROJECT_PATH}/${sample}/outs/possorted_genome_bam.bam" "${LOCAL_DIR}/${BAM_NAME}"
    cp "${PROJECT_PATH}/${sample}/outs/possorted_genome_bam.bam.bai" "${LOCAL_DIR}/${BAM_NAME}.bai"
    BAM_FOUND=true
fi

if [ "$BAM_FOUND" = false ]; then
    echo "Error: BAM file for sample ${sample} not found in expected locations." >&2
    exit 1
fi

cd "${LOCAL_DIR}"
echo "Successfully staged BAM file."

# Index the BAM file if the index doesn't exist
if [ ! -f "${BAM_NAME}.bai" ]; then
    echo "BAM index not found. Indexing BAM file..."
    samtools index -@ "${SLURM_CPUS_PER_TASK}" "${BAM_NAME}"
fi

# 2. Create BAM List
echo "--- 2. Creating BAM list file ---"
echo "${project}_${sample},${BAM_NAME}" > bam.lst

# 3. Monopogen Preprocessing
echo "--- 3. Running Monopogen PreProcess ---"
SAMP_DIR="${project}_${sample}"
mkdir -p "${SAMP_DIR}"

python "${MONOPOGEN_PATH}/src/Monopogen.py" preProcess \
    -b bam.lst \
    -t "${SLURM_CPUS_PER_TASK}" \
    -a "${MONOPOGEN_PATH}/apps" \
    -o "${SAMP_DIR}"

echo "Preprocessing complete."

# 4. Monopogen Germline Variant Calling
echo "--- 4. Running Monopogen Germline Calling ---"
# Note: Monopogen's germline step might use more cores than allocated by -c. 
# Here we use 8 as specified in the original script. Adjust if needed.
python "${MONOPOGEN_PATH}/src/Monopogen.py" germline \
    -t 8 \
    -a "${MONOPOGEN_PATH}/apps" \
    -r "${MONOPOGEN_PATH}/resource/GRCh38.region.50MB.lst" \
    -o "${SAMP_DIR}" \
    -p "${IMPUTATION_PATH}/" \
    -g "${REFERENCE_PATH}" \
    -m 3 \
    -s all

echo "Germline calling complete."

# 5. Copy Results Back
echo "--- 5. Copying results back to project directory ---"
MONOPOGEN_OUTPUT_DIR="${PROJECT_PATH}/monopogen"
mkdir -p "${MONOPOGEN_OUTPUT_DIR}"
cp -r "${SAMP_DIR}" "${MONOPOGEN_OUTPUT_DIR}/"

echo "Results successfully copied to ${MONOPOGEN_OUTPUT_DIR}"
echo "--- Monopogen run finished successfully! ---"

# The 'trap cleanup EXIT' will handle the removal of the local directory.
