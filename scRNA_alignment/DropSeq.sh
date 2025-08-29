#!/bin/bash
# ==============================================================================
# Shell Script for DropSeq Single-Cell RNA-Seq Analysis
# ==============================================================================


source ~/.bashrc

tissue=$1
project=$2

dropseq_dir="/path/to/Drop-seq_tools-2.5.4"
cd /path/to/scRNA/${tissue}/raw_data/${project}

module load jdk/11.0.10
module load picard
module load star

#-------------------------------------
# Make reference for drop-seq alignment
#-------------------------------------
#bash ${dropseq_dir}/create_Drop-seq_reference_metadata.sh  \
#	-n dropseq-ref \
#	-s human \
#	-g /path/to/RhapRef/gencodev29-20181205.gtf \
#        -d ${dropseq_dir} \
#	-r /path/to/genome_reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#	-o /path/to/DropseqRef \
#	-a /path/to/STAR \
#	-v 

sample=`cat sample.list | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1`
F1=${sample}_S1_L001_R1_001.fastq.gz
F2=${sample}_S1_L001_R2_001.fastq.gz

#-------------------------------------
# convert fastq to  unaligned bam
#-------------------------------------
picard FastqToSam F1=$F1 F2=$F2 \
	OUTPUT=${sample}_unaligned.bam \
	SAMPLE_NAME=${sample}

#-------------------------------------
# QC, trimming, and STAR alignment
#-------------------------------------
mkdir -p ${sample}_alignment

bash ${dropseq_dir}/Drop-seq_alignment.sh\
	-g /path/to/DropseqRef/STAR \
	-r /path/to/DropseqRef/dropseq-ref.fasta.gz \
	-d ${dropseq_dir} \
	-s /path/to/star-2.7.9a/bin/Linux_x86_64_static/STAR \
	-o ${sample}_alignment \
	-t ${sample}_alignment \
	${sample}_unaligned.bam

#-------------------------------------
# Call gene count matrix from aligned bam
#-------------------------------------
if [ ! -d raw_matrix ]; then
	mkdir raw_matrix
fi

if [ -f ${sample}_alignment/final.bam ]; then 
	${dropseq_dir}/DigitalExpression I=${sample}_alignment/final.bam O=raw_matrix/${project}_${sample}.dge.txt.gz SUMMARY=raw_matrix/${project}_${sample}.dge.summary.txt NUM_CORE_BARCODES=100
fi
#
