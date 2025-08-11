#!/bin/bash
##############################################################################
# Functional annotation using snpEff for all SNPs tested in sc-eQTL analysis
##############################################################################
#SBATCH --mem=10G 
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -J snpeff
#SBATCH --output=LOG/snpeff_%A.out
#SBATCH --error=LOG/snpeff_%A.err

module load bedtools
module load jdk
snpEff=/path/to/snpEff/snpEff.jar
snpSift=/path/to/snpEff/SnpSift.jar
vcfPL=/path/to/snpEff/scripts/vcfEffOnePerLine.pl


tissue=$1
ancestry=$2
sig_path="/path/to/scRNA/${tissue}/sc-eQTL/results/sig"
geno_path="/path/to/scRNA/${tissue}/sc-eQTL/genotype"
func_anno_path="/path/to/scRNA/${tissue}/sc-eQTL/results/functional_annotation/snpEff"

if [ ! -d $func_anno_path ]; then
  mkdir -p $func_anno_path
fi

dir_local=/data/$$
mkdir -p $dir_local

cp ${sig_path}/*_${ancestry}_pc5_esnp_5en8.bed $dir_local
cp ${geno_path}/${tissue}_${ancestry}_indivpruned.bim $dir_local 

cd $dir_local 
mkdir -p overlaps

################# Tissue - Ancestry level #####################
### create vcf for null SNP
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${tissue}_${ancestry}_null.vcf
awk '{ 
    # Extract columns from BIM file
    chrom=$1; id=$2; pos=$4; ref=$6; alt=$5; 
    # Write to VCF format
    printf "%s\t%s\t%s\t%s\t%s\t.\tPASS\t.\n", chrom, pos, id, ref, alt; 
}' ${tissue}_${ancestry}_indivpruned.bim >> ${tissue}_${ancestry}_null.vcf 

rm *bed

### Running snpEff
# null
java -Xmx8g -jar $snpEff -v GRCh38.86 ${tissue}_${ancestry}_null.vcf > ${tissue}_${ancestry}_null_snpeff.vcf
cat ${tissue}_${ancestry}_null_snpeff.vcf | $vcfPL | java -jar $snpSift extractFields - ID 'ANN[*].EFFECT' | tail -n +2 | sort | uniq > ${tissue}_${ancestry}_null_snpeff.txt

# Overlapping of annotations in eSNPs
# Overlapping of annotations in null
cat ${tissue}_${ancestry}_null_snpeff.txt | cut -f2 | sort | uniq -c > all_anno_null.txt
for anno in `cat all_anno_null.txt | awk '$1>1' | awk '{print $2}'`; do
  awk -v n=$anno '$2==n' ${tissue}_${ancestry}_null_snpeff.txt  | cut -f 1 > overlaps/${tissue}_${ancestry}_null_${anno}_overlaps.txt 
  echo "$anno: $(wc -l overlaps/${tissue}_${ancestry}_null_${anno}_overlaps.txt  | cut -d' ' -f1) overlapping eSNPs"
done

# Copy back
cp overlaps/* $func_anno_path
cp *snpeff.vcf $func_anno_path

cd /data
rm -rf $dir_local
