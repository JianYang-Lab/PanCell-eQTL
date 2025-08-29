#!/bin/bash
# ==============================================================================
# Shell Script for demutiplexing using demuxlet or freemuxlet
# ==============================================================================


set -e

# tissue, project, SLURM array task ID, and method are passed as arguments
tissue="$1"
project="$2"
i="$3"
demulti_method="$4" # demuxlet or freemuxlet

# Define base directory for the project
dir="/path/to/scRNA/${tissue}/raw_data/${project}"

# Get sample name and expected number of individuals (N) from a metadata file
# Using SLURM_ARRAY_TASK_ID ($i) to select the correct line
sample=$(cat "${dir}/sample.n" | head -n "${i}" | tail -n 1 | awk '{print $1}')
N=$(cat "${dir}/sample.n" | head -n "${i}" | tail -n 1 | awk '{print $2}')

# Create a unique temporary directory on fast local storage
dir_local="/data/$$"
mkdir -p "$dir_local"
cd "$dir_local"

# Ensure cleanup happens even if the script fails
trap 'cd /data; rm -rf "$dir_local"' EXIT

if [ "$demulti_method" == "freemuxlet" ]; then
    echo "Running freemuxlet for sample ${sample}"
    # Copy required files to the local directory
    cp "${dir}/${sample}/outs/possorted_genome_bam.bam"* .
    cp "${dir}/${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" .
    cp /path/to/resource/1000genomes/1KGP3_3202_GRCh38_filtered_phased_sorted.vcf.gz .

    # Define file variables for clarity
    bam_file="possorted_genome_bam.bam"
    dbsnp_vcf="1KGP3_3202_GRCh38_filtered_phased_sorted.vcf.gz"
    group_list="barcodes.tsv.gz"

    # 1. Pile-up against a common SNP VCF
    popscle dsc-pileup --sam "${bam_file}" --vcf "${dbsnp_vcf}" --group-list "${group_list}" --out pileup

    # 2. Demultiplex without genotypes, providing the expected number of samples
    popscle freemuxlet --plp pileup --out freemuxlet_result --nsample "${N}"

    # Copy results back, excluding large input files
    for f in $(ls | grep -v bam | grep -v barcodes); do
        cp -r "$f" "${dir}/${sample}/"
    done

elif [ "$demulti_method" == "demuxlet" ]; then 
    echo "Running demuxlet for sample ${sample}"
    # Copy required files to the local directory
    cp "${dir}/${sample}/outs/possorted_genome_bam.bam"* .
    cp "${dir}/${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" .
    cp /path/to/merged.vcf . # VCF file for your specific multiplexed individuals

    # Define file variables
    bam_file="possorted_genome_bam.bam"
    group_list="barcodes.tsv.gz"
    known_vcf="merged.vcf"

    # 1. Pile-up against the VCF of the KNOWN multiplexed individuals
    # Using known_vcf instead of a generic population VCF
    popscle dsc-pileup --sam "${bam_file}" --vcf "${known_vcf}" --group-list "${group_list}" --out pileup

    # 2. Demultiplex using the pileup and known genotypes
    popscle demuxlet --plp pileup --vcf "${known_vcf}" --field GT --out demuxlet_result

    # Copy results back
    for f in $(ls | grep -v bam | grep -v barcodes | grep -v merged.vcf); do
        cp -r "$f" "${dir}/${sample}/"
    done
fi
