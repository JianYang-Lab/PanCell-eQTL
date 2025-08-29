#!/bin/bash
#------------------------------------------------------------
# dropEst Alignment Pipeline Script
# Uses STAR for alignment and dropEst for expression estimation
# Suitable for platforms including:
# 10x cel_seq2 config_desc drop_seq_velocyto drop_seq iclip 
# indrop_v1_2 indrop_v3 seq_well split_seq 
#------------------------------------------------------------

# Exit immediately if a command exits with a non-zero status
set -e

#------------------------------------------------------------
# Part 0: Data preparation
#------------------------------------------------------------
module load tophat
module load samtools

tissue=$1
project=$2
i=$3
config=$4 
# 10x cel_seq2 config_desc drop_seq_velocyto drop_seq iclip indrop_v1_2 indrop_v3 seq_well split_seq
dir=/path/to/scRNA/${tissue}/raw_data/${project}
sample=`cat ${dir}/sample.list | head -n ${i} | tail -n 1`
dir_local=/data/$$

mkdir -p $dir_local

cp ${dir}/${sample}_S*.fastq.gz $dir_local

cd $dir_local
#------------------------------------------------------------
# Part 1: Set Up Variables and Paths
#------------------------------------------------------------

# User-configurable variables

# Paths to input FASTQ files (gzip compressed)
READ1_FASTQ=`ls ${sample}*R1_001.fastq.gz`
READ2_FASTQ=`ls ${sample}*R2_001.fastq.gz`

# Reference genome files
GENOME_DIR=/path/to/genome_reference/STAR/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v36_oh100
GENOME_FASTA=/path/to/genome_reference/GRCh38_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa
ANNOTATION_GTF=/path/to/genome_reference/GRCh38_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf

# Output directories
ALIGNMENT_DIR=./
DROP_EST_OUTPUT_DIR=${sample}

# Number of threads to use
THREADS=8

# Sample name
SAMPLE_NAME=${sample}

# STAR parameters (advanced users may customize)
MAX_MEMORY="160G"          # Adjust based on your system's RAM

# dropEst parameters (advanced users may customize)
DROPEST_BIN=/path/to/dropest  # Path to dropEst executable
DROPTAG_BIN=/path/to/dropEst/build/droptag
DROPEST_CONFIG="/path/to/dropEst/configs/${config}.xml"

#------------------------------------------------------------
# Part 2: dropTag fastq for alignment
#------------------------------------------------------------
if [ ! -f ${dir}/bam/${project}_${sample}.bam ]; then
$DROPTAG_BIN -c $DROPEST_CONFIG $READ1_FASTQ $READ2_FASTQ -p 8
cat ${READ2_FASTQ}.tagged.*.fastq.gz > ${sample}_tagged.fastq.gz

#------------------------------------------------------------
# Part 3: Perform Alignment with STAR
#------------------------------------------------------------

echo "Starting alignment with STAR..."

# Create output directory for alignment
#mkdir -p "$ALIGNMENT_DIR"
BowtieIndex=/path/to/Bowtie_index/GRCh38_noalt_as/GRCh38_noalt_as
# Run TopHat aligner
tophat2 -p $THREADS \
       	--no-coverage-search -g 1 \
	-G ${ANNOTATION_GTF} \
	-o ${ALIGNMENT_DIR} \
       	${BowtieIndex} ${sample}_tagged.fastq.gz

samtools index accepted_hits.bam -@ 8

cp accepted_hits.bam ${dir}/bam/${project}_${sample}.bam
cp accepted_hits.bam.bai ${dir}/bam/${project}_${sample}.bam.bai
tail align_summary.txt -n 1
echo "Alignment completed. BAM file generated at $ALIGNMENT_DIR"

# Define the path to the aligned BAM file
ALIGNED_BAM="$ALIGNMENT_DIR/accepted_hits.bam"

else
# Define the path to the aligned BAM file
ALIGNED_BAM="${dir}/bam/${project}_${sample}.bam"
fi

#------------------------------------------------------------
# Part 4: Process BAM with dropEst
#------------------------------------------------------------

echo "Starting dropEst processing..."

# Create output directory for dropEst
mkdir -p "$DROP_EST_OUTPUT_DIR"

# dropEst configuration file for inDrop (provide your own or use default)
DROPEST_CONFIG="/path/to/dropEst/configs/${config}.xml"  # Adjust path as needed

# Run dropEst
$DROPEST_BIN -V -m -g $ANNOTATION_GTF -c $DROPEST_CONFIG \
    -o ${DROP_EST_OUTPUT_DIR}/${SAMPLE_NAME}.dropEst \
    $ALIGNED_BAM

cp -r $DROP_EST_OUTPUT_DIR $dir

cd /data
rm -rf $dir_local

echo "dropEst processing completed. Output generated at $DROP_EST_OUTPUT_DIR"

#------------------------------------------------------------
# Completion Message
#------------------------------------------------------------

echo "Pipeline execution finished successfully!"
