#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perform permutations for type I error checking.           
"""

import os
import sys
from datetime import datetime
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import numpy as np
from collections import defaultdict
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas {pd.__version__}")

def load_genotype(tissue,ancestry):
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
    #other_cols = ['Age_pred', 'Gender_pred']
    covariates_df = covariates_df[pf_cols+pc_cols+cd_cols+other_cols] 
    covariates_df = covariates_df.loc[phenotype_df.columns]
    return(phenotype_df, phenotype_pos_df, covariates_df)

def permute_by_batch(phenotype_df):
     # Extract batch information from column names
    pro_batch = [col.split("_")[0] for col in phenotype_df.columns]
    # Create batch to columns mapping
    batch_to_columns = defaultdict(list)
    for col, batch in zip(phenotype_df.columns, pro_batch):
        batch_to_columns[batch].append(col)
    # Make a copy of the dataframe for permutation
    permuted_df = phenotype_df.copy()
    # Identify batches with single samples 
    single_sample_batches = [b for b, cols in batch_to_columns.items() if len(cols) == 1]
    multi_sample_batches = [b for b, cols in batch_to_columns.items() if len(cols) > 1]
    for batch in multi_sample_batches:
        cols = batch_to_columns[batch]
        n = len(cols)
        # Generate permutation (avoid identity)
        while True:
            perm = np.random.permutation(n)
            if not np.array_equal(perm, np.arange(n)):
                break
        # Get permuted columns and rename to original
        permuted_cols = [cols[i] for i in perm]
        rename_dict = {pc: oc for pc, oc in zip(permuted_cols, cols)}
        batch_data = phenotype_df[permuted_cols].rename(columns=rename_dict)
        # Update columns in permuted_df
        permuted_df[cols] = batch_data[cols]
    # Process pooled single-sample batches
    pooled_cols = [col for b in single_sample_batches for col in batch_to_columns[b]]
    if len(pooled_cols) > 1:
        # Generate permutation (avoid identity)
        n = len(pooled_cols)
        while True:
            perm = np.random.permutation(n)
            if not np.array_equal(perm, np.arange(n)):
                break 
        # Get permuted columns and rename to original
        permuted_pooled = [pooled_cols[i] for i in perm]
        rename_dict = {pc: oc for pc, oc in zip(permuted_pooled, pooled_cols)}
        pooled_data = phenotype_df[permuted_pooled].rename(columns=rename_dict)
        # Update columns in permuted_df
        permuted_df[pooled_cols] = pooled_data[pooled_cols]
    # Verification
    changes = (phenotype_df != permuted_df).sum().sum()
    total = phenotype_df.size
    print(f"Changed {changes}/{total} values ({100*changes/total:.2f}%)")
    return(permuted_df)
    
 
if __name__=='__main__':
    # Read in cell type and ancestry
    #ct_list = ['CD4TNC','CD4TEM','Treg','CD8TNC','CD8TEM','CD8TEMRA','MAIT','BIN','BMem','DC','MonoC','MonoNC','NKn','NKp','pDC','Plasma']
    tissue = sys.argv[1]
    ancestry = sys.argv[2]
    celltype = sys.argv[3]
    pc_num = int(sys.argv[4])
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Ancestry:{ancestry} Cell type:{celltype}  nPC:{pc_num}', flush=True)
    
    os.chdir(f"/path/to/{tissue}/sc-eQTL")

    # load genotypes and variants into data frames
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load genotype for ancestry {ancestry} started', flush=True)
    genotype_df,variant_df = load_genotype(tissue, ancestry)
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load genotype for ancestry {ancestry} completed', flush=True)

    # load genes.gtf gene_name, gene_chr, gene_start, gene_end, strand
    gene_gtf = load_genes_gtf()
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load genes.gtf completed', flush=True)

    for i in [2]:
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] permutation {i} starts', flush=True)
        prefix = f"results/perm/{celltype}_{ancestry}_pc{pc_num}_perm{i}" 
        if not os.path.isfile(f"{prefix}.cis_qtl_pairs.1.parquet"): 
            # load phenotypes and covariates
            phenotype_df, phenotype_pos_df, covariates_df = load_phenotype_and_covariates(ancestry, celltype, pc_num)
            #permuted_df = permute_by_batch(phenotype_df)
            # Shuffle phenotype columns
            shuffled_columns = np.random.permutation(phenotype_df.columns)
            phenotype_df.columns = shuffled_columns
            covariates_df = covariates_df.loc[shuffled_columns]
            print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Load phenotypes for {ancestry}_{celltype}_{pc_num} and shuffle y completed', flush=True)
            
            # cis-QTL mapping: summary statistics for all variant-phenotype pairs
            cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                            prefix, covariates_df, output_dir='.', maf_threshold=0.01)

            print('['+datetime.now().strftime("%b %d %H:%M:%S")+f'] Completed eQTL mapping for {ancestry}_{celltype}_pc{pc_num}_perm{i}', flush=True) 
                 
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Session completed', flush=True) 
