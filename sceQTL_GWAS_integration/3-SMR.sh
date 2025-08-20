#!/bin/bash
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -c 8
#SBATCH --mem=36G
#SBATCH --qos=huge
#SBATCH --job-name=sceQTL_GWAS_SMR
#SBATCH --output=sceQTL_GWAS_SMR_%j.out
#SBATCH --error=sceQTL_GWAS_SMR_%j.err

# SMR mode: default, diff, fdr, multi
# Usage: sbatch 3-SMR.sh <mode> <GWAS_id> <GWAS_ancestry> <tissue> <ancestry> <celltype>
# Example
# sbatch 3-SMR.sh default 1 EAS Blood EAS B
# sbatch 3-SMR.sh diff 1 EAS Blood EUR CD4TNC

mode=$1 

# GWAS info
GWAS_id=$2
GWAS_ancestry=$3

# eQTL info
tissue=$4
ancestry=$5
celltype=$6

export GWAS_id
export tissue
export ancestry
export celltype
export GWAS_ancestry

# 1. Load paths ---------------------------------------------------
# SMR software path
SMR=/path/to/smr-1.3.1-linux-x86_64/smr
# Ref path
ref_dir=/path/to/resource/1000genomes
# GWAS path
if [ $GWAS_ancestry == "EAS" ]; then
	GWAS_list=/path/to/GWAS_${GWAS_ancestry}_list.txt
else
	GWAS_list=/path/to/GWAS_${GWAS_ancestry}_list.txt
fi
GWAS_trait=`head -n ${GWAS_id} $GWAS_list | tail -n 1 | cut -f 1`
GWAS_path=`head -n ${GWAS_id} $GWAS_list | tail -n 1 | cut -f 4`
GWAS_path_local=`basename $GWAS_path`
# eQTL path
eQTL_path=/path/to/${tissue}/sc-eQTL/results/besd/${tissue}_${celltype}_${ancestry}_sceQTL
snp_probe_path=/path/to/${tissue}/sc-eQTL/results/sig/${tissue}_${celltype}_${ancestry}_snp_probe_fdr05.list
# Out path
if [ $mode == "default" ]; then
    out_dir=/path/to/${tissue}/sceqtl_gwas_integration/SMR/raw
elif [ $mode == "diff" ]; then
    out_dir=/path/to/${tissue}/sceqtl_gwas_integration/SMR/raw_diff
elif [ $mode == "fdr" ]; then
    out_dir=/path/to/${tissue}/sceqtl_gwas_integration/SMR/raw_fdr
elif [ $mode == "multi" ]; then
    out_dir=/path/to/${tissue}/sceqtl_gwas_integration/SMR/raw_multi
fi
if [ ! -d ${out_dir} ]; then
    mkdir -p $out_dir
fi


# 2. cp data to /data ----------------------------------------------
dir_local=/data/$$
mkdir -p $dir_local

cp ${GWAS_path} $dir_local
cp ${eQTL_path}* $dir_local

cd $dir_local

# 3. run SMR analysis for each trait and eQTL ----------------------------------------------
if [[ $GWAS_path_local == *.gz ]]; then
    gunzip $GWAS_path_local
    GWAS_summary=`echo $GWAS_path_local | sed 's/\.gz//'`
else 
    GWAS_summary="$GWAS_path_local"
fi

echo -e "[Start SMR] GWAS: ${GWAS_trait} sceQTL: ${tissue}_${celltype}_${ancestry}"
for ((i=1;i<=22;i++)); do
 #i=2
    echo -e "Processing chr${i}"
    if [ ! -f ${out_dir}/SMR_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i}.smr ]; then
    ref=${ref_dir}/${ancestry}/1KGP_${ancestry}_chr${i}_maf01_mind05_hwe1en6_msite05_hg38
    cp ${ref}* $dir_local
    if [ $mode == "default" ]; then
        ${SMR} --bfile 1KGP_${ancestry}_chr${i}_maf01_mind05_hwe1en6_msite05_hg38 \
            --gwas-summary ${GWAS_summary} \
            --beqtl-summary ${tissue}_${celltype}_${ancestry}_sceQTL \
            --maf 0.01 \
            --thread-num 8 \
            --out SMR_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i}
    elif [ $mode == "diff" ]; then
	 ${SMR} --bfile 1KGP_${ancestry}_chr${i}_maf01_mind05_hwe1en6_msite05_hg38 \
            --gwas-summary ${GWAS_summary} \
            --beqtl-summary ${tissue}_${celltype}_${ancestry}_sceQTL \
            --maf 0.01 \
            --thread-num 8 \
            --out SMR_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i} \
	    --diff-freq-prop 0.5
    elif [ $mode == "fdr" ]; then
        ${SMR} --bfile 1KGP_${ancestry}_chr${i}_maf01_mind05_hwe1en6_msite05_hg38 \
            --extract-target-snp-probe ${snp_probe_path} \
            --gwas-summary ${GWAS_summary} \
            --beqtl-summary ${tissue}_${celltype}_${ancestry}_sceQTL \
            --maf 0.01 \
            --peqtl-smr 0.05 \
            --thread-num 8 \
            --out SMR_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i}
    elif [ $mode == "multi" ]; then
        ${SMR} --bfile 1KGP_${ancestry}_chr${i}_maf01_mind05_hwe1en6_msite05_hg38 \
            --gwas-summary ${GWAS_summary} \
            --beqtl-summary ${tissue}_${celltype}_${ancestry}_sceQTL \
            --maf 0.01 \
            --smr-multi  \
            --thread-num 8 \
            --out SMR_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i}
    fi 
    cp SMR_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i}* $out_dir
    rm SMR_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i}* 
    fi
done

cd /data
rm -rf $dir_local

if [ $? -ne 0 ]; then
    echo "Error raised"
    exit 1
fi

echo "All jobs completed successfully for SMR on GWAS: ${GWAS_trait} sceQTL: ${tissue}_${celltype}_${ancestry}"
