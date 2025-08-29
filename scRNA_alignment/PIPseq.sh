#!/bin/bash
# ==============================================================================
# Shell Script for PipSeq Single-Cell RNA-Seq Analysis
# ==============================================================================


tissue=$1
project=$2
i=$3

cd /path/to/scRNA/${tissue}/raw_data/${project}

echo $project

SRR=`cat sample.list | head -n $i | tail -n 1`
mkdir $SRR

pipseeker_path="/path/to/pipseeker"
star_index_path="/path/to/genome_reference/STAR/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v36_oh100"

${pipseeker_path}/pipseeker full --chemistry v5  --fastq="${SRR}" \
	--sorted-bam \
	--star-index-path ${star_index_path} \
	--annotation ${pipseeker_path}/human-pbmc-references/references/human-pbmc-v4.csv \
	--output-path ${SRR} \
	--threads 8
