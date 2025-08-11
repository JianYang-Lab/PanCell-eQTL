#!/bin/bash
##############################################################################
# Description: Functional annotation using chromHMM for all SNPs tested in 
# sc-eQTL analysis.
##############################################################################
#SBATCH --mem=10G 
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -J bedtools
#SBATCH --output=LOG/bedtools_%A.out
#SBATCH --error=LOG/bedtools_%A.err

module load bedtools

tissue=$1
ancestry=$2
sig_path="/path/to/${tissue}/sc-eQTL/results/sig"
geno_path="/path/to/${tissue}/sc-eQTL/genotype"
chromHMM_path="/storage/yangjianLab/chenchang/resource/chromHMM/"
chromHMM_dir=$3
func_anno_path="/path/to/${tissue}/sc-eQTL/results/functional_annotation/chromHMM"

dir_local=/data/$$
mkdir -p $dir_local

cp -r ${chromHMM_path}/${chromHMM_dir} $dir_local
cp ${geno_path}/${tissue}_${ancestry}_indivpruned.bim $dir_local 

cd $dir_local 
mkdir -p overlaps

################# Tissue - Ancestry level #####################
### all SNPs tested (Null)
awk 'BEGIN {FS="\t"; OFS="\t"} \
  {chrom = "chr"$1; \
  start = $4 - 1; \
  print chrom, start, $4, $2 \
  }' ${tissue}_${ancestry}_indivpruned.bim > ${tissue}_${ancestry}_null.bed

for annotation in ${chromHMM_dir}/*.bed; do
  base=$(basename $annotation .bed)
  bedtools intersect -a ${tissue}_${ancestry}_null.bed -b $annotation -wa -u > overlaps/${tissue}_${ancestry}_null_${chromHMM_dir}_${base}_overlaps.bed
  echo "$base: $(wc -l overlaps/${tissue}_${ancestry}_null_${chromHMM_dir}_${base}_overlaps.bed | cut -d' ' -f1) overlapping eSNPs"
done

cp overlaps/*bed $func_anno_path

cd /data
rm -rf $dir_local
