#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@file: 4-tensorqtl.py
@desc: eQTL mapping using tensorQTL
"""

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
    plink_prefix_path = f"genotype/{tissue}_{ancestry}_indivpruned"
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
    return(phenotype_df, phenotype_pos_df, covariates_df)

def process_parquet_chr(celltype, ancestry, pc_num, chr, egene_ids, threshold_dict):
    parquet_file = f'results/raw/{celltype}_{ancestry}_pc{pc_num}.cis_qtl_pairs.{chr}.parquet'
    df = pd.read_parquet(parquet_file)
    df = df[df['phenotype_id'].isin(egene_ids)]
    df = df[(df['af'] >= 0.01) & (df['af'] <= 0.99)] #! filter out MAF<0.01
    sig_indice = df['pval_nominal'] < df['phenotype_id'].apply(lambda x: threshold_dict[x])
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

# select and output significant eQTLs
def select_sig_eqtl(celltype, ancestry, pc_num, fdr = 0.05):
    # permutation result
    permutation_tsv = f"results/raw/{celltype}_{ancestry}_pc{pc_num}_perm.tsv"
    gene_df = pd.read_csv(permutation_tsv, sep='\t', index_col=0)
    #gene_info = pd.merge(gene_gtf, gene_df, left_index=True, right_index=True, how='right') 
    # eGenes (apply FDR threshold)
    if sum(gene_df.columns.str.contains('pval_nominal_threshold'))>0:
        egene_df = gene_df.loc[gene_df['qval']<=fdr, ['pval_nominal_threshold', 'pval_nominal', 'pval_beta']].copy() 
        egene_df.rename(columns={'pval_nominal': 'min_pval_nominal'}, inplace=True)
        egene_ids = set(egene_df.index)
        threshold_dict = egene_df['pval_nominal_threshold'].to_dict() 
        # load parquet files
        signif_df = []
        for i in range(1,23):
            signif_df.append(process_parquet_chr(celltype, ancestry, pc_num, i, egene_ids, threshold_dict))
        signif_df = pd.concat(signif_df, axis=0)
        signif_df = signif_df.merge(egene_df, left_on='phenotype_id', right_index=True)
        out_path = f"results/sig/{celltype}_{ancestry}_pc{pc_num}_sig.tsv.gz"
        signif_df.to_csv(out_path, sep='\t', index=False, float_format='%.6g', compression='gzip')
        n_egenes = signif_df.phenotype_id.nunique()
        n_eqtls = signif_df.shape[0] 
    else: 
        n_genes = 0
        n_eqtls = 0
    record_file = 'results/sig/eQTL_summary.tsv'
    append_eqtl_summary(record_file, celltype, ancestry, pc_num, n_egenes, n_eqtls)

if __name__=='__main__':
    # Read in tissue, cell type and ancestry
    if len(sys.argv) != 5:
        print("Usage: python 4-tensorqtl.py <tissue> <ancestry> <celltype> <pc_num>")
        sys.exit(1)
        
    tissue = sys.argv[1]
    ancestry = sys.argv[2]
    celltype = sys.argv[3]
    pc_num = int(sys.argv[4])

    # tissue is set in the bash script
    os.chdir(f"/path/to/{tissue}/sc-eQTL")

    # load genotypes and variants into data frames
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load genotype for ancestry {ancestry} started', flush=True)
    genotype_df,variant_df = load_genotype(ancestry)
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load genotype for ancestry {ancestry} completed', flush=True)

    # load genes.gtf gene_name, gene_chr, gene_start, gene_end, strand
    #gene_gtf = load_genes_gtf()
    #print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load genes.gtf completed', flush=True)

    #for celltype in ct_list:
    prefix = f"results/raw/{celltype}_{ancestry}_pc{pc_num}" 
    if not os.path.isfile(f"{prefix}_perm.tsv"): 
        # load phenotypes and covariates
        phenotype_df, phenotype_pos_df, covariates_df = load_phenotype_and_covariates(ancestry, celltype, pc_num)
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load phenotypes for {ancestry}_{celltype}_{pc_num} completed', flush=True) 
        
        # cis-QTL mapping: permutations: controlling for global FDR < 0.05
        cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, maf_threshold=0.01)
        tensorqtl.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)
        cis_df.to_csv(f"{prefix}_perm.tsv", sep='\t')
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Completed permutations for {ancestry}_{celltype}_{pc_num}', flush=True)
       
        # Conditional independent analysis 
        indep_df = cis.map_independent(genotype_df, variant_df, cis_df,
                                        phenotype_df, phenotype_pos_df, covariates_df, maf_threshold=0.01)
        indep_df.to_csv(f"{prefix}_indep.tsv", sep='\t', index=False)
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Completed conditional eQTL mapping for {ancestry}_{celltype}_{pc_num}', flush=True)
             
        # cis-QTL mapping: summary statistics for all variant-phenotype pairs
        cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                        prefix, covariates_df, output_dir='.', maf_threshold=0.01)

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Completed eQTL mapping for {ancestry}_{celltype}_{pc_num}', flush=True)
        
    if os.path.isfile(f"{prefix}_perm.tsv"):
        if not os.path.isfile(f"results/sig/{celltype}_{ancestry}_pc{pc_num}_sig.tsv.gz"):
            # Select and output significant eQTLs
            select_sig_eqtl(celltype, ancestry, pc_num, fdr = 0.05)
            print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Completed annotation', flush=True)
           
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Session completed', flush=True) 
