#!/bin/bash
#SBATCH --mem=20G
#SBATCH --qos=normal
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH -c 1
#SBATCH -J LDSC

# ==============================================================================
# LDSC-SEG Pipeline for sc-eQTL and GWAS
# ==============================================================================

# 1. Define Paths and Variables
# ------------------------------------------------------------------------------
TISSUE=$1
CT=$2
GWAS_NAME=$3

WORKDIR="/path/to/your/workdir"
EQTL_DIR="${WORKDIR}/scRNA/${TISSUE}/sc-eQTL_R/results/query"
LDSC_RES="/path/to/LDSC_resource"

# Set output directories
OUT_DIR="${WORKDIR}/scRNA/${TISSUE}/sceqtl_gwas_integration_R/LDSC_SEG_results/"
MUNGED_GWAS="/path/to/GWAS/munged/${GWAS_NAME}_EUR_munged"
mkdir -p ${OUT_DIR}

# P-value threshold for defining eGenes (Standard genome-wide significance)
PVAL_THRESH=5e-8

echo "=========================================================="
echo " Step 2-4: Process Cell-Type Specific eQTLs"
echo "=========================================================="

echo "----------------------------------------------------------"
echo " Processing Cell Type: ${CT}"
echo "----------------------------------------------------------"
    
EQTL_FILE="${EQTL_DIR}/${TISSUE}_${CT}_EUR_sceQTL_query.txt"
GENESET="${OUT_DIR}/annotations_fdr05/${CT}.GeneSet"
    
# --------------------------------------------------------------------------
# Step 2: Extract Significant eGenes
# --------------------------------------------------------------------------
# Columns: 1:SNP 2:Chr 3:BP 4:A1 5:A2 6:Freq 7:Probe 8:Probe_Chr 9:Probe_bp 
#          10:Gene 11:Orientation 12:b 13:se 14:p
if [ ! -f ${GENESET} ]; then
    awk -v thresh=${PVAL_THRESH} 'NR>1 && $14 < thresh {print $7}' ${EQTL_FILE} | sort | uniq > ${GENESET}
fi
GENE_COUNT=$(wc -l < ${GENESET})
echo "Extracted ${GENE_COUNT} significant eGenes for ${CT}."
    
# --------------------------------------------------------------------------
# Step 3: Make Annotations and Compute LD Scores (Loop over chromosomes)
# --------------------------------------------------------------------------
# Standard SEG window is 100kb around the gene boundary
JOB_IDS=()
for CHR in {1..22}; do
    ANNOT_PREFIX="${OUT_DIR}/annotations_fdr05/${CT}.${CHR}"
    
    # Check if files already exist to skip
    if [[ -f "${ANNOT_PREFIX}.l2.ldscore.gz" ]]; then continue; fi

    # Create a small temporary shell script for this specific chromosome
    JOB_CMD="python /path/to/ldsc/make_annot.py \
        --gene-set-file ${GENESET} \
        --gene-coord-file /path/to/resource/ENSG_coord.txt \
        --windowsize 100000 \
        --bimfile ${LDSC_RES}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR}.bim \
        --annot-file ${ANNOT_PREFIX}.annot.gz; \ 
    python /path/to//ldsc/ldsc.py \
        --thin-annot --l2 \
        --print-snps ${LDSC_RES}/hapmap3_snps/w_hm3.snpslist \
        --bfile ${LDSC_RES}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR} \
        --ld-wind-cm 1 \
        --annot ${ANNOT_PREFIX}.annot.gz \
        --out ${ANNOT_PREFIX}"

    JID=$(qsubshcom "$JOB_CMD" 1 10G LDSC_score_${TISSUE}_${CT}_${CHR} 100000:00:00 "--qos=huge -p intel-sc3,amd-ep2,amd-ep2-short")
    JOB_IDS+=($JID)
done

# CRITICAL: Wait for all qsubshcom jobs submitted by this session to finish
# Note: Ensure your local environment's 'wait' or a specific check loop is used.
# Since qsubshcom usually returns immediately, we check for file completion:

echo "Waiting for all chromosome jobs to complete..."
if [ ${#JOB_IDS[@]} -gt 0 ]; then
    echo "Waiting for background jobs to finish: ${JOB_IDS[*]}"
    
    while true; do
        STILL_RUNNING=0
        for ID in "${JOB_IDS[@]}"; do
            # Use qstat or squeue depending on your scheduler. 
            # If squeue (Slurm):
            if squeue -j $ID 2>/dev/null | grep -q $ID; then
                STILL_RUNNING=$((STILL_RUNNING + 1))
            fi
        done
        
        if [ "$STILL_RUNNING" -eq 0 ]; then
            break
        fi
        
        echo "Still waiting for $STILL_RUNNING jobs..."
        sleep 120
    done
fi
echo "LD Scores computed for ${CT}."

# --------------------------------------------------------------------------
# Step 4: Partition Heritability
# --------------------------------------------------------------------------
# We test your custom annotation *jointly* with the baseline v2.2 model 
# to control for general genomic architectures (like coding regions, enhancers, etc).
# The output will tell us if your eGenes explain MORE heritability than expected.

echo "Partitioning Heritability for ${CT}..."
python /path/to//ldsc/ldsc.py \
    --h2 ${MUNGED_GWAS}.sumstats.gz \
    --ref-ld-chr ${OUT_DIR}/annotations_fdr05/${CT}.,${LDSC_RES}/baselineLD_v2.2/baselineLD. \
    --w-ld-chr ${LDSC_RES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --frqfile-chr ${LDSC_RES}/1000G_Phase3_frq/1000G.EUR.QC. \
    --overlap-annot \
    --out ${OUT_DIR}/res_${GWAS_NAME}_${CT}

echo "Completed analysis for ${CT}."
