#!/bin/bash

# ==============================================================================
# Shell Script for BD Rhapsody Single-Cell RNA-Seq Analysis
#
# This script runs the BD Rhapsody WTA (Whole Transcriptome Analysis) pipeline.
# It takes paired-end FASTQ files as input, performs read trimming, aligns
# reads to a reference genome, corrects UMIs, and generates a gene
# expression matrix.
#
# ==============================================================================

# --- Exit immediately if a command exits with a non-zero status ---
set -e
# --- Treat unset variables as an error ---
set -u
# --- Ensure that pipes report the first error ---
set -o pipefail

# --- Usage Message ---
usage() {
    echo "Usage: $0 -1 <READ1_FASTQ> -2 <READ2_FASTQ> -s <SAMPLE_NAME> -o <OUTPUT_DIR>"
    echo "  -1  Path to Read 1 FASTQ file (contains cell barcodes and UMIs)."
    echo "  -2  Path to Read 2 FASTQ file (contains transcript sequence)."
    echo "  -s  A unique name for the sample (e.g., 'Sample_A_T-cells')."
    echo "  -o  Path to the directory where all output will be saved."
    exit 1
}

# --- Parse Command-Line Arguments ---
while getopts "1:2:s:o:" opt; do
    case ${opt} in
        1 ) r1_fastq=$OPTARG;;
        2 ) r2_fastq=$OPTARG;;
        s ) sample_name=$OPTARG;;
        o ) output_dir=$OPTARG;;
        \? ) usage;;
    esac
done

# Check if all required arguments are provided
if [ -z "${r1_fastq-}" ] || [ -z "${r2_fastq-}" ] || [ -z "${sample_name-}" ] || [ -z "${output_dir-}" ]; then
    echo "ERROR: Missing required arguments."
    usage
fi

# ==============================================================================
# --- CONFIGURATION ---
#! IMPORTANT: Update these paths to match your system's environment.
# ==============================================================================

# --- Path to the main BD Rhapsody pipeline executable ---
# This might be a shell script or a direct call to a Java/Python program.
BD_PIPELINE_EXECUTABLE="/path/to/your/bd-rhapsody-pipeline/main.sh"

# --- Reference Genome and Annotation Files ---
# The reference must be indexed with STAR.
REFERENCE_DIR="/path/to/your/reference/GRCh38_BD_STAR_v1.0"
REFERENCE_GTF="/path/to/your/reference/Homo_sapiens.GRCh38.108.gtf"

# --- (Optional) Supplemental Reference Files ---
# Use these if you have AbSeq (antibody) or Sample Tag data.
# Leave the string empty ("") if not in use.
ABSEQ_REFERENCE_CSV="" # e.g., "/path/to/abseq_reference.csv"
SAMPLE_TAGS_REFERENCE_CSV="" # e.g., "/path/to/sample_tags.csv"

# --- Resource Allocation ---
THREADS=16
MEMORY_LIMIT_GB=120 # Max memory for the pipeline in Gigabytes

# ==============================================================================
# --- SCRIPT EXECUTION ---
# ==============================================================================

echo "🚀 Starting BD Rhapsody Pipeline for sample: $sample_name"
echo "--------------------------------------------------------"
echo "Read 1: $r1_fastq"
echo "Read 2: $r2_fastq"
echo "Output Directory: $output_dir"
echo "Reference GTF: $REFERENCE_GTF"
echo "--------------------------------------------------------"
start_time=$(date +%s)

# --- 1. Setup Output Directory ---
echo "⚙️  Creating output directory..."
mkdir -p "$output_dir"

# --- 2. Construct the Main Pipeline Command ---
# The core command to run the BD analysis.
# We build it in an array to cleanly handle optional arguments.
CMD=(
    "$BD_PIPELINE_EXECUTABLE"
    --r1 "$r1_fastq"
    --r2 "$r2_fastq"
    --sample-name "$sample_name"
    --output-dir "$output_dir"
    --reference-path "$REFERENCE_DIR"
    --annotation-path "$REFERENCE_GTF"
    --exact-cell-count 0 # Set to a number to keep only the top N cells, or 0 to auto-detect
    --threads $THREADS
    --max-memory-gb $MEMORY_LIMIT_GB
)

# Add optional arguments if the files are specified
if [ -n "$ABSEQ_REFERENCE_CSV" ]; then
    echo "INFO: AbSeq reference provided."
    CMD+=("--abseq-reference" "$ABSEQ_REFERENCE_CSV")
fi

if [ -n "$SAMPLE_TAGS_REFERENCE_CSV" ]; then
    echo "INFO: Sample Tags reference provided for multiplexing."
    CMD+=("--sample-tags-reference" "$SAMPLE_TAGS_REFERENCE_CSV")
fi

# --- 3. Execute the Pipeline ---
echo "🔬 Running Analysis. This may take a long time..."
echo
# Print the full command that will be executed for reproducibility
echo "Full command:"
echo "${CMD[@]}"
echo

# Execute the command
"${CMD[@]}"

end_time=$(date +%s)
runtime=$((end_time - start_time))

# --- 4. Finalization ---
echo
echo "--------------------------------------------------------"
echo "✅ Pipeline finished successfully!"
echo "Total runtime: $((runtime / 3600))h $(((runtime / 60) % 60))m $((runtime % 60))s"
echo "Find your results in: $output_dir"
echo "--------------------------------------------------------"