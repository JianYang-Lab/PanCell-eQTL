#!/bin/bash
#SBATCH --mem=1G
#SBATCH --qos=huge
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --job-name=sc-eQTL_pipeline


# Check user inputs

if [ $# -ne 3 ]; then

    echo "Usage: $0 <tissue> <cell_type> <ancestry>"
    exit 1

fi

ts=$1
ct=$2
anc=$3

# ----------- Job Submission Functions -----------
submit_normalization() {
    local output_file="/path/to/celltype_matrix/norm_rint_tsv/${ts}_expr_norm_rint_${ct}_${anc}.tsv"
    
    if [ ! -f "$output_file" ]; then

        local input_file="../h5ad_celltype/${ts}_rawQC_${ct}_${anc}.h5ad"
        local target_size=$(du -BG "$input_file" | cut -d G -f 1)
        local apply_size=$((target_size * 10))
        
        # Determine resources

        if [ $apply_size -lt 450 ]; then

            local partition="intel-sc3,amd-ep2,amd-ep2-short"
            local qos="normal"
        else

            local partition="intel-fat"
            local qos="hmem"
        fi

        
        # Submit normalization job with environment activation

        sbatch --job-name=norm_${ct}_${anc} \
               --output=LOG/norm_${ct}_${anc}.out \
               --error=LOG/norm_${ct}_${anc}.err \
               --mem=${apply_size}G \
               --partition=$partition \
               --qos=$qos \
               --parsable \
               --wrap "eval \"\$(micromamba shell hook -s bash)\" && micromamba run -n scanpy_env python 1-normalize_rint.py $ts $anc $ct"
    fi

}

submit_peer() {
    local peer_file="/path/to/covariates/peer_factors/${ct}_${anc}_peer_factors.tsv"
    
    #if [ -f "/path/to/celltype_matrix/norm_rint_tsv/${ts}_expr_norm_rint_${anc}_${ct}.tsv" ] && [ ! -f "$peer_file" ]; then
    if [ ! -f "$peer_file" ]; then
        local input_file="/path/to/celltype_matrix/norm_rint_tsv/${ts}_expr_norm_rint_${anc}_${ct}.tsv"
        local target_size=$(du -BG "$input_file" | cut -d G -f 1)
        local apply_size=$((target_size * 4))

        # Determine resources

        if [ $apply_size -lt 450 ]; then
            local partition="intel-sc3,amd-ep2"
            local qos="huge"
        else
            local partition="intel-fat"
            local qos="hmem"
        fi

        
        local dependency_option=""
        if [ -n "$1" ]; then

            dependency_option="--dependency=afterok:$1"
        fi

        # Submit PEER job with environment activation

        sbatch --job-name=peer_${ct}_${anc} \
               --output=LOG/peer_${ct}_${anc}.out \
               --error=LOG/peer_${ct}_${anc}.err \
               --mem=${apply_size}G \
               --partition=$partition \
               --qos=$qos \
               $dependency_option \
               --parsable \
               --wrap "/path/to/anaconda3/envs/peer/bin/Rscript 2-peer.R $ts $ct $anc"
    fi

}

submit_prep() {
    local output_file="/path/to/covariates/final/${ct}_${anc}_cov_final.txt"
    local output_file1="/path/to/celltype_matrix/bed/${ct}_${anc}_expr.bed"

    if [ ! -f "$output_file" ] || [ ! -f "${output_file1}" ]; then

        local dependency_option=""
        if [ -n "$1" ]; then

            dependency_option="--dependency=afterok:$1"
        fi
 
        # Submit data prep job with environment activation

        sbatch --job-name=prep_${ct}_${anc} \
               --output=LOG/prep_${ct}_${anc}.out \
               --error=LOG/prep_${ct}_${anc}.err \
               --qos=huge \
	       --mem=10G \
               --partition=intel-sc3,amd-ep2,amd-ep2-short \
               $dependency_option \
               --parsable \
               --wrap "eval \"\$(micromamba shell hook -s bash)\" && micromamba run -n R_env Rscript 3-prepare_data_for_eqtl.R $ts $ct $anc"
    fi

}

submit_tensorqtl() {
    # Submit tensorQTL job with environment activation

    local dependency_option=""
    if [ -n "$1" ]; then

        dependency_option="--dependency=afterok:$1"
    fi

    sbatch --job-name=tensorqtl_${ct}_${anc} \
           --output=LOG/tensorqtl_${ct}_${anc}.out \
           --error=LOG/tensorqtl_${ct}_${anc}.err \
           --partition=a40-tmp,a40-quad,v100 \
	   --qos=gpu-huge \
	   --gres=gpu:1 \
           --mem=60G \
           $dependency_option \
           --wrap "eval \"\$(micromamba shell hook -s bash)\" && micromamba run -n tensorqtl python 4-tensorqtl.py $ts $anc $ct 5"
}

# ----------- Main Pipeline -----------
echo "Starting pipeline for $ct $anc"
cd /path/to/${ts}/sc-eQTL/

module load R/4.3.1
# 1. Submit normalization
norm_jid=$(submit_normalization)
echo "Submitted normalization job: $norm_jid"

# Wait for normalization to complete and check for output file
if [ -n "$norm_jid" ]; then
    echo "Waiting for normalization job $norm_jid to complete..."
    while squeue -j $norm_jid > /dev/null 2>&1; do
        sleep 60 # Wait for 1 minute before checking again
    done
    echo "Normalization job $norm_jid completed."
fi

# 2. Submit PEER
peer_jid=$(submit_peer $norm_jid)
echo "Submitted PEER job: $peer_jid"

# Wait for PEER to complete and check for output file
if [ -n "$peer_jid" ]; then
    echo "Waiting for PEER job $peer_jid to complete..."
    while squeue -j $peer_jid > /dev/null 2>&1; do
        sleep 60 # Wait for 1 minute before checking again
    done
    echo "PEER job $peer_jid completed."
fi

# 3. Submit data prep
prep_jid=$(submit_prep $peer_jid)
echo "Submitted data prep job: $prep_jid"

# Wait for data prep to complete and check for output file
if [ -n "$prep_jid" ]; then
    echo "Waiting for data prep job $prep_jid to complete..."
    while squeue -j $prep_jid > /dev/null 2>&1; do
        sleep 60 # Wait for 1 minute before checking again
    done
    echo "Data prep job $prep_jid completed."
fi

# 4. Submit tensorQTL
tensor_jid=$(submit_tensorqtl $prep_jid)
echo "Submitted tensorQTL job: $tensor_jid"

# Wait for tensorQTL to complete
if [ -n "$tensor_jid" ]; then
    echo "Waiting for tensorQTL job $tensor_jid to complete..."
    while squeue -j $tensor_jid > /dev/null 2>&1; do
        sleep 60 # Wait for 1 minute before checking again
    done
    echo "TensorQTL job $tensor_jid completed."
fi

echo "Pipeline completed for $ts $ct $anc."
