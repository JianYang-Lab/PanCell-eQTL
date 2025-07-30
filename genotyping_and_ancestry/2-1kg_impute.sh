#!/bin/bash
#SBATCH --mem=2G
#SBATCH --qos=normal
#SBATCH -c 2
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --job-name=submit_imputation
#SBATCH --output=imputation_submitter_%j.log

# ==============================================================================
# Two-step Imputation Submission Script: Part I
# ==============================================================================
# This script submits per-sample genotype imputation jobs to a cluster
# using Minimac4 and a 1000 Genomes reference panel.
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error.
set -o pipefail # Exit status of a pipeline is the last command to exit with a non-zero status.

# --- Function Definitions ---

# Prints usage information for the script.
usage() {
    echo "Usage: $0 --tissue <tissue> --project <project_id> --sample <sample_name>"
    echo ""
    echo "Options:"
    echo "  --tissue    The name of the tissue (e.g., Colon)."
    echo "  --project   The ID of the project to process (e.g., PRJNA477814)."
    echo "  --sample    The name of the sample to process (e.g., Sample1)."
    echo "  -h, --help  Display this help message."
    exit 1
}

# --- Configuration & Argument Parsing ---

# Static paths for software and reference data
readonly MINIMAC4_EXEC="/storage/yangjianLab/chenchang/SOFTWARE/Minimac4/bin/minimac4"
readonly REF_PATH="/storage/yangjianLab/chenchang/resource/1000genomes/minimac_ref"
readonly OUTPUT_DIR="../VCF_1KGimputed/VCF_persample"

# Parse command-line arguments
if [ "$#" -eq 0 ]; then
    usage
fi

while [ "$#" -gt 0 ]; do
    case "$1" in
        --tissue)   tissue="$2"; shift 2;;
        --project)  project_id="$2"; shift 2;;
        --sample)   sample_name="$2"; shift 2;;
        -h|--help)  usage;;
        *)          echo "Unknown option: $1" >&2; usage;;
    esac
done

# Check if all required arguments were provided
if [ -z "${tissue:-}" ] || [ -z "${project_id:-}" ] || [ -z "${sample_name:-}" ]; then
    echo "Error: Missing one or more required arguments." >&2
    usage
fi

# --- Environment Setup ---
# Ensure necessary environment variables or modules are loaded.
# 'source ~/.bashrc' is general; be more specific if possible.
source ~/.bashrc

# --- Main Workflow ---

echo "--- Starting Imputation Job Submission for Project: ${project_id}, Sample: ${sample_name} ---"

# Construct the base directory using the 'tissue' argument
readonly BASE_DIR="/path/to/your/project/${tissue}/raw_data"

# Create the main output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

readonly project_dir="${BASE_DIR}/${project_id}"
if [ ! -d "${project_dir}/monopogen" ]; then
    echo "Error: Monopogen directory not found for project ${project_id} at ${project_dir}/monopogen" >&2
    exit 1
fi

# Find the specific sample directory
sample_dir="${project_dir}/monopogen/${sample_name}"
if [ ! -d "${sample_dir}" ]; then
    echo "Error: Sample directory ${sample_name} not found in ${project_dir}/monopogen" >&2
    exit 1
fi

vcf_path="${sample_dir}/vcf"
echo "Processing sample: ${sample_name}"

if [ ! -d "${vcf_path}" ]; then
    echo "Warning: VCF directory not found for sample ${sample_name}. Skipping."
    exit 0 # Exit gracefully as there's nothing to do for this sample
fi

# Iterate through chromosomes 1 to 22
for i in {1..22}; do
    phased_vcf="${vcf_path}/chr${i}.phased.sorted.vcf.gz"
    imputed_vcf="${OUTPUT_DIR}/${sample_name}_chr${i}_imputed.vcf.gz"
    ref_msav="${REF_PATH}/1KGP3_GRCh38_chr${i}.msav"
    job_name="IPT_${project_id}_${sample_name}_${i}"

    # Check if the output file already exists to avoid re-running completed jobs
    if [ -f "${imputed_vcf}" ]; then
        echo "  Output for chr${i} already exists. Skipping."
        continue
    fi

    # Check if required input files exist
    if [ ! -f "${phased_vcf}" ]; then
        echo "  Warning: Phased VCF for chr${i} not found. Skipping."
        continue
    fi
    if [ ! -f "${ref_msav}" ]; then
        echo "  Error: Reference .msav file for chr${i} not found. Cannot proceed." >&2
        exit 1
    fi

    # Define the command to be submitted
    impute_command="$MINIMAC4_EXEC ${ref_msav} ${phased_vcf} -o ${imputed_vcf} --format GT --threads=2 --min-ratio-behavior skip"

    # Submit the job using the qsubshcom wrapper
    echo "  Submitting job for chr${i}..."
    qsubshcom "${impute_command}" 2 8G "${job_name}" "100:00:00" "-p intel-sc3,amd-ep2,amd-ep2-short --qos=huge"
done


echo "--- All jobs for sample ${sample_name} in project ${project_id} have been submitted. ---"
