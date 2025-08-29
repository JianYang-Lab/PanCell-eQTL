#!/bin/bash

# ==============================================================================
# Smart-seq Analysis Pipeline adapted from HCA
#
# This script performs a complete analysis for a single Smart-seq sample,
# including:
#   1. Alignment of FASTQ reads to a reference genome using HISAT2.
#   2. Quality control of the resulting alignment using Picard Tools.
#   3. Gene expression quantification using RSEM on the aligned BAM file.
#
# The script is designed to be portable and robust. All external file paths
# are defined in the configuration section, and inputs are passed via
# command-line arguments.
#
# ==============================================================================


# --- Stop on any error ---
set -e
set -o pipefail

# --- Usage Message ---
usage() {
    echo "Usage: $0 -p <PROJECT_DIR> -s <SAMPLE_PREFIX> -e <PAIRED_END_FLAG>"
    echo "  -p  Path to the project directory containing FASTQ files."
    echo "  -s  The sample prefix (e.g., 'ERR3648608.lite.2')."
    echo "  -e  Flag indicating paired-end data ('true' or 'false')."
    exit 1
}

# --- Parse Command-Line Arguments ---
while getopts "p:s:e:" opt; do
    case ${opt} in
        p ) project_dir=$OPTARG;;
        s ) sample_prefix=$OPTARG;;
        e ) paired_end=$OPTARG;;
        \? ) usage;;
    esac
done

# Check if all required arguments are provided
if [ -z "$project_dir" ] || [ -z "$sample_prefix" ] || [ -z "$paired_end" ]; then
    usage
fi

# ==============================================================================
# --- CONFIGURATION ---
# IMPORTANT: All paths to reference files and software should be defined here.
# ==============================================================================

# --- Software Dependencies (for documentation purposes) ---
# This script requires:
# - hisat2 (~2.2.1)
# - samtools (~1.15)
# - picard (can be a .jar file) (~2.27)
# - rsem (~1.3.3)
# - Java (~8 or higher)

# --- Reference Genome and Annotation Files ---
GENOME_REF_FASTA="/path/to/genome_reference/GRCh38/GRCh38.fa"
RRNA_INTERVALS="/path/to/genecode/gencode.v44.rRNA.interval_list"
GENE_REF_FLAT="/path/to/genecode/refFlat.txt"

# --- Index Files ---
HISAT2_REF_INDEX="/path/to/hisat2_reference/GRCh38_hisat2/GRCh38_hisat2"
RSEM_REF_INDEX="/path/to/rsem_reference/rsem_hisat_trans_ref_index"

# --- Tool Paths ---
PICARD_JAR="/path/to/picard.jar"

# --- Library strandedness ---
# Options: NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND
STRANDEDNESS="NONE"

# --- Resources ---
THREADS=8


# ==============================================================================
# --- SETUP ---
# ==============================================================================

echo "Starting analysis for sample: $sample_prefix"
echo "Project directory: $project_dir"
echo "Paired-end: $paired_end"
echo "----------------------------------------"

cd "$project_dir"

# Define input FASTQ files
if [ "$paired_end" = "true" ]; then
    fastq1="${sample_prefix}_1.fastq.gz"
    fastq2="${sample_prefix}_2.fastq.gz"
    # Check if files exist
    if [ ! -f "$fastq1" ] || [ ! -f "$fastq2" ]; then
        echo "ERROR: Paired-end FASTQ files not found for prefix $sample_prefix"
        exit 1
    fi
else
    fastq1="${sample_prefix}.fastq.gz"
    if [ ! -f "$fastq1" ]; then
        echo "ERROR: Single-end FASTQ file not found: $fastq1"
        exit 1
    fi
fi

# Define output names
output_prefix="${sample_prefix}_out"
qc_prefix="${output_prefix}_qc"
rsem_prefix="${output_prefix}_rsem"

# Define final BAM filename
sorted_bam="${qc_prefix}.bam"


# ==============================================================================
# --- STEP 1: ALIGNMENT (HISAT2) ---
# Align FASTQ reads to the reference genome.
# The output is piped directly to samtools to create a sorted, indexed BAM file.
# This avoids writing a large intermediate unsorted BAM file.
# ==============================================================================

echo "STEP 1: Aligning with HISAT2..."

# Prepare common HISAT2 options
hisat2_opts=(
    -p $THREADS
    -x "$HISAT2_REF_INDEX"
    --seed 12345
    --rg-id "$sample_prefix" --rg "SM:$sample_prefix" --rg "LB:$sample_prefix"
    --rg "PL:ILLUMINA" --rg "PU:$sample_prefix"
    --new-summary --summary-file "${qc_prefix}.hisat2.summary.log"
    --met-file "${qc_prefix}.hisat2.met.txt" --met 5
    -k 10
    --secondary
)

# Run alignment and pipe to sorting/indexing
if [ "$paired_end" = "true" ]; then
    hisat2 "${hisat2_opts[@]}" -1 "$fastq1" -2 "$fastq2" | \
    samtools sort -@ $THREADS -O bam -o "$sorted_bam" -
else
    # Note: These parameters are highly specific and disable spliced alignment.
    # This may be desired for certain analyses but is not typical for standard RNA-seq.
    hisat2 "${hisat2_opts[@]}" \
        -U "$fastq1" \
        --mp 1,1 \
        --np 1 \
        --score-min L,0,-0.1 \
        --no-mixed \
        --no-softclip \
        --no-discordant \
        --rdg 99999999,99999999 \
        --rfg 99999999,99999999 \
        --no-spliced-alignment | \
    samtools sort -@ $THREADS -O bam -o "$sorted_bam" -
fi

echo "Indexing BAM file..."
samtools index "$sorted_bam"

echo "HISAT2 alignment complete. Sorted BAM created: $sorted_bam"
echo "----------------------------------------"


# ==============================================================================
# --- STEP 2: QUALITY CONTROL (Picard) ---
# Collect a comprehensive set of QC metrics from the aligned BAM file.
# ==============================================================================
echo "STEP 2: Running Picard for Quality Control..."

# Collect multiple metrics with a single command
java -jar "$PICARD_JAR" CollectMultipleMetrics \
    VALIDATION_STRINGENCY=SILENT \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    I="$sorted_bam" \
    O="${qc_prefix}" \
    FILE_EXTENSION=".txt" \
    PROGRAM=null \
    PROGRAM=CollectAlignmentSummaryMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    PROGRAM=CollectGcBiasMetrics \
    PROGRAM=CollectBaseDistributionByCycle \
    PROGRAM=QualityScoreDistribution \
    PROGRAM=MeanQualityByCycle \
    PROGRAM=CollectSequencingArtifactMetrics \
    PROGRAM=CollectQualityYieldMetrics \
    REFERENCE_SEQUENCE="$GENOME_REF_FASTA" \
    ASSUME_SORTED=true

# Collect specific RNA-seq metrics
java -Xmx5g -jar "$PICARD_JAR" CollectRnaSeqMetrics \
    VALIDATION_STRINGENCY=SILENT \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    INPUT="$sorted_bam" \
    OUTPUT="${qc_prefix}.rna_metrics.txt" \
    REF_FLAT="$GENE_REF_FLAT" \
    RIBOSOMAL_INTERVALS="$RRNA_INTERVALS" \
    STRAND_SPECIFICITY="$STRANDEDNESS" \
    CHART_OUTPUT="${qc_prefix}.rna.coverage.pdf"

echo "Picard QC complete."
echo "----------------------------------------"

# --- Optional: Mark Duplicates ---
# This step is commented out by default but can be enabled if needed.
# It identifies and flags or removes PCR duplicates.
# Removing duplicates is sometimes debated for RNA-seq, as it can remove
# legitimate reads from highly expressed genes.
#
# echo "STEP 2b: Marking duplicates..."
# java -Xmx5g -XX:ParallelGCThreads=$THREADS -jar "$PICARD_JAR" MarkDuplicates \
#     VALIDATION_STRINGENCY=SILENT \
#     INPUT="$sorted_bam" \
#     OUTPUT="${qc_prefix}.dedup.bam" \
#     ASSUME_SORTED=true \
#     METRICS_FILE="${qc_prefix}.duplicate_metrics.txt" \
#     REMOVE_DUPLICATES=false # Set to true to remove duplicates instead of just flagging


# ==============================================================================
# --- STEP 3: QUANTIFICATION (RSEM) ---
# Calculate gene and isoform expression counts.
# CRITICAL FIX: This now uses the BAM file generated in Step 1,
# avoiding a redundant second alignment.
# ==============================================================================
echo "STEP 3: Quantifying expression with RSEM..."

rsem_opts=(
    --bam # Use the BAM file as input
    -p $THREADS
    --calc-pme
    --seed 555
    --single-cell-prior
    --time
)

if [ "$paired_end" = "true" ]; then
    rsem-calculate-expression --paired-end "${rsem_opts[@]}" "$sorted_bam" "$RSEM_REF_INDEX" "$rsem_prefix"
else
    rsem-calculate-expression "${rsem_opts[@]}" "$sorted_bam" "$RSEM_REF_INDEX" "$rsem_prefix"
fi

echo "RSEM quantification complete."
echo "----------------------------------------"


# ==============================================================================
# --- STEP 4: CLEANUP & ORGANIZE ---
# Move all output files into a structured directory for the sample.
# ==============================================================================
echo "STEP 4: Organizing output files..."

mkdir -p "${output_prefix}/QC"

# Move QC metrics and plots
mv ${qc_prefix}*txt "${output_prefix}/QC/"
mv ${qc_prefix}*pdf "${output_prefix}/QC/"

# Move logs and main BAM file
mv ${qc_prefix}.hisat2.summary.log "${output_prefix}/"
mv ${sorted_bam}* "${output_prefix}/" # Moves .bam and .bai

# Move RSEM results
mv ${rsem_prefix}.* "${output_prefix}/"

echo "Analysis complete for $sample_prefix!"
echo "All results are in the '${output_prefix}' directory."
# ==============================================================================
