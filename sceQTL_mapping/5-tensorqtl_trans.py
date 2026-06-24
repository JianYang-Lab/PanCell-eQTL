##############################################################################
# Script information
# Step 5: TensorQTL
# Author: Chang Chen
# Date: 2024-09-26
# Description: Perform sc-eQTL mapping using TensorQTL
##############################################################################
import os
import sys
from datetime import datetime
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas {pd.__version__}")

def load_genotype(ancestry):
    plink_prefix_path = f"genotype/Blood_{ancestry}_indivpruned_allcissig_mappability_filtered"
    pr = genotypeio.PlinkReader(plink_prefix_path)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    return(genotype_df, variant_df)
    
def load_genes_gtf(annotation_gtf = "/path/to/resource/genes.gtf"): 
    gene_gtf = pd.read_csv(annotation_gtf, sep='\t', index_col=1)
    gene_gtf.columns = ['gene_name','gene_chr','gene_start','gene_end','strand'] 
    return(gene_gtf)
    
def load_phenotype_and_covariates(ancestry, celltype, pc_num):
    expression_bed = f"celltype_matrix/bed/{celltype}_{ancestry}_expr.bed"
    covariates_file = f"covariates/final/{celltype}_{ancestry}_cov_final.txt"
    # Load phenotype
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
    pf_cols = list(covariates_df.columns[covariates_df.columns.str.startswith('PF')])
    pc_cols = [f'PC{i}' for i in range(1, pc_num + 1)]
    cd_cols = list(covariates_df.columns[covariates_df.columns.str.startswith('Condition_group')])
    other_cols = ['Age_pred', 'Gender_pred']
    covariates_df = covariates_df[pf_cols+pc_cols+cd_cols+other_cols] 
    covariates_df = covariates_df.loc[phenotype_df.columns]
    # ! Filter out low mappability genes
    gene_mappability = pd.read_csv("/path/to/resource/Mappability/hg38_gene_mappability.txt.gz",sep='\t',header=None)
    gene_mappability.columns = ["phenotype_id","mappability"]
    gene_mappability['phenotype_id'] = [ph.split(".")[0] for ph in gene_mappability['phenotype_id']]
    gene_keep = gene_mappability.loc[gene_mappability['mappability']>=0.8,"phenotype_id"]
    phenotype_df = phenotype_df.loc[phenotype_df.index.isin(gene_keep)] #! filter out genes with mappability<0.8
    phenotype_pos_df = phenotype_pos_df.loc[phenotype_pos_df.index.isin(gene_keep)] #! filter out genes with mappability<0.8 
    return(phenotype_df, phenotype_pos_df, covariates_df)

def process_parquet_chr(celltype, ancestry, pc_num, chr, egene_ids, threshold_dict):
    parquet_file = f'results/raw/{celltype}_{ancestry}_pc{pc_num}.cis_qtl_pairs.{chr}.parquet'
    df = pd.read_parquet(parquet_file)
    df = df[df['phenotype_id'].isin(egene_ids)]
    sig_indice = df['pval_nominal']<df['phenotype_id'].apply(lambda x: threshold_dict[x])
    signif_df = df[sig_indice]
    print(f'Chr{chr} processed') 
    return(signif_df)

def append_eqtl_summary(record_file, celltype, ancestry, pc_num, n_egenes, n_eqtls):
    file_exists = os.path.isfile(record_file)
    header = "cell_type\tancestry\tnPC\tneGenes\tneQTLs\n"
    data_row = f"{celltype}\t{ancestry}\t{pc_num}\t{n_egenes}\t{n_eqtls}\n"
    with open(record_file, 'a' if file_exists else 'w') as file:
        if not file_exists:
            file.write(header)
        file.write(data_row)

if __name__=='__main__':
    # Read in cell type and ancestry
    ancestry = sys.argv[1]
    celltype = sys.argv[2]
    pc_num = int(sys.argv[3])
    os.chdir("/path/to/Blood/sc-eQTL_R")

    # load genotypes and variants into data frames
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load genotype for ancestry {ancestry} started', flush=True)
    genotype_df,variant_df = load_genotype(ancestry)
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load genotype for ancestry {ancestry} completed', flush=True)

    prefix = f"results/raw/{celltype}_{ancestry}_pc{pc_num}" 
    # load phenotypes and covariates
    phenotype_df, phenotype_pos_df, covariates_df = load_phenotype_and_covariates(ancestry, celltype, pc_num)
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load phenotypes for {ancestry}_{celltype}_{pc_num} completed', flush=True) 
    
    if not os.path.isfile(f"{prefix}_allcissig_trans_r2.parquet"): 
            # trans-eQTL mapping
            trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, batch_size=10000,
                           return_sparse=True, pval_threshold=0.05, maf_threshold=0.01)
            # remove cis-associations
            trans_df = trans.filter_cis(trans_df, phenotype_pos_df, variant_df, window=5000000)
            trans_df.to_parquet(f"{prefix}_allcissig_trans_r2.parquet")
            print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Completed trans eQTL mapping for {ancestry}_{celltype}_{pc_num}', flush=True)
    else:
            trans_df = pd.read_parquet(f"{prefix}_allcissig_trans_r2.parquet")
   
    # 1. Run permutations (perm_series is a pd.Series of global parameters)
    perm_series = trans.map_permutations(genotype_df, covariates_df, batch_size=10000, maf_threshold=0.01)

    # 2. Apply permutations to top gene-snp pairs
    top_df = trans_df.loc[trans_df.groupby('phenotype_id')['pval'].idxmin()].copy()
    trans.apply_permutations(perm_series, top_df) 
   
    # 3. Calculate Q-value 
    tensorqtl.calculate_qvalues(top_df, fdr=0.05, qvalue_lambda=0.85) 
    top_df.to_csv(f"{prefix}_trans_perm.tsv", sep="\t")
    
    egenes_df = top_df[top_df['qval'] < 0.05]
    
    if len(egenes_df) > 0:
        global_nominal_threshold = egenes_df['pval_nominal_threshold'].max()
        significant_trans_df = trans_df[(trans_df["phenotype_id"].isin(egenes_df["phenotype_id"])) & (trans_df["pval"]<=global_nominal_threshold)] 
        significant_trans_df.to_parquet(f"{prefix}_significant_trans_eqtls.parquet")
    else:
        print("No significant trans-eGenes found at FDR < 0.05.")
        
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Session completed', flush=True) 
