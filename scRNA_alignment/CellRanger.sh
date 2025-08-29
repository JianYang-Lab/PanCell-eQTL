#!/bin/bash
# ==============================================================================
# Shell Script for 10X Single-Cell RNA-Seq Analysis
# ==============================================================================

module load cellranger

ref="/path/to/genome_reference/GRCh38_10X/refdata-gex-GRCh38-2020-A"
tissue=$1
project=$2
i=$3
CR_dir=/path/to/scRNA/${tissue}/raw_data/${project}

sample=`cat ${CR_dir}/sample.list | head -n $i | tail -n 1`
SRR=`ls ${CR_dir} | grep fastq.gz | rev | cut -d _ -f 5- | rev | sort | uniq | grep $sample | paste -sd ',' -`

dir_local=/data/$$
mkdir -p $dir_local
		
cp ${CR_dir}/${sample}*fastq.gz $dir_local
		
cd $dir_local	

cellranger count --id=${sample} \
		--create-bam=true \
		--transcriptome=${ref} \
		--fastqs=./ \
        --sample=${SRR} \
		--localcores=12 \
        --localmem=56

if [ $? -ne 0 ]; then
	echo "CellRanger failed"
	exit 1
fi

cp -r ${sample} $CR_dir
cd /data

rm -rf $dir_local


