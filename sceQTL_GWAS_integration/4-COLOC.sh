#!/bin/bash
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH --qos=huge
#SBATCH --job-name=sceQTL_COLOC
#SBATCH --output=sceQTL_COLOC_%j.out
#SBATCH --error=sceQTL_COLOC_%j.err

module load R/4.3.1

# COLOC mode: default, diff, or fdr
# default: use ancestry-matched GWAS summary statistics
# diff: use ancestry-mismatched GWAS summary statistics
# fdr: use FDR-corrected eQTL probes
mode=$1

# GWAS info
GWAS_id=$2
GWAS_ancestry=$3

# eQTL info
tissue=$4
ancestry=$5
celltype=$6

# 1. Load paths -----------------------------------------------------
SMR=/path/to/smr-1.3.1-linux-x86_64/smr
COLOC=COLOC_perm.R
GWAS_list=/path/to/GWAS_${GWAS_ancestry}_list.txt
GWAS_trait=`head -n ${GWAS_id} $GWAS_list | tail -n 1 | cut -f 1`
GWAS_path=`head -n ${GWAS_id} $GWAS_list | tail -n 1 | cut -f 4`
GWAS_path_local=`basename $GWAS_path`
GWAS_N=`head -n ${GWAS_id} $GWAS_list | tail -n 1 | cut -f 5`
if [ -z "$GWAS_N" ]; then
	GWAS_N="NA"
fi
eQTL_path=/path/to/${tissue}/sc-eQTL/results/besd/${tissue}_${celltype}_${ancestry}_sceQTL
QTLN=`cat /path/to/crossTissue_analysis/tissue_celltype_ancestry_ss_new.txt | grep $tissue | grep $ancestry | awk -v ct=$celltype '$3==ct' | awk '{print $4}'`
bed_dir="/path/to/${tissue}/sc-eQTL/celltype_matrix/bed"
probe_path=/path/to/${tissue}/sc-eQTL/results/sig

# Out path
if [ $mode == "default" ]; then
    out_dir=/path/to/${tissue}/sceqtl_gwas_integration/COLOC/raw
elif [ $mode == "diff" ]; then
    out_dir=/path/to/${tissue}/sceqtl_gwas_integration/COLOC/raw_diff
elif [ $mode == "fdr" ]; then
    out_dir=/path/to/${tissue}/sceqtl_gwas_integration/COLOC/raw_fdr
fi
if [ ! -d ${out_dir} ]; then
    mkdir -p $out_dir
fi

# 2. cp data to /data ----------------------------------------------
dir_local=/data/$$
mkdir -p $dir_local

cp ${GWAS_path} $dir_local
cp ${eQTL_path}* $dir_local
cp ${bed_dir}/${celltype}_${ancestry}_expr.bed $dir_local

cd $dir_local

# 3. run COLOC analysis for each trait and eQTL --------------------------------------------
echo -e "[Start COLOC] GWAS: ${GWAS_trait}_${GWAS_ancestry} sceQTL: ${tissue}_${celltype}_${ancestry}"
for ((i=1;i<=22;i++)); do
if [ ! -f ${out_dir}/COLOC_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i}.coloc ]; then
	
	if [ "$mode" == "default" ] || [ "$mode" == "diff" ]; then
		Rscript ${COLOC} \
			--SMR $SMR \
			--sumstats ${GWAS_path_local} \
			--GWASN ${GWAS_N} \
			--besd ${tissue}_${celltype}_${ancestry}_sceQTL \
			--QTLName ${tissue}_${celltype}_${ancestry} \
			--QTLN ${QTLN} \
			--perm FALSE \
			--expr_bed ${celltype}_${ancestry}_expr.bed \
			--chr ${i} \
			--out ${out_dir}/COLOC_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i} \
			--tmp_dirt $dir_local	
	elif [ $mode == "fdr" ]; then	
		cp ${probe_path}/${tissue}_${celltype}_${ancestry}_snp_probe_fdr05.list $dir_local
		Rscript ${COLOC} \
			--SMR $SMR \
			--sumstats ${GWAS_path_local} \
			--GWASN ${GWAS_N} \
			--besd ${tissue}_${celltype}_${ancestry}_sceQTL \
			--QTLName ${tissue}_${celltype}_${ancestry} \
			--QTLN ${QTLN} \
			--expr_bed ${celltype}_${ancestry}_expr.bed \
			--perm TRUE \
			--probe_path ${tissue}_${celltype}_${ancestry}_snp_probe_fdr05.list \
			--chr ${i} \
			--out ${out_dir}/COLOC_${GWAS_trait}_${tissue}_${celltype}_${ancestry}_${GWAS_ancestry}_chr${i} \
			--tmp_dirt $dir_local
	fi
fi
done

cd /data
rm -rf $dir_local

if [ $? -ne 0 ]; then
    echo "Error raised"
    exit 1
fi

echo "All jobs completed successfully for COLOC on GWAS: ${GWAS_trait} sceQTL: ${celltype}_${ancestry}"

