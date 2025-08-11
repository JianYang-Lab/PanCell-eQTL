#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script reads in single-cell expression data, normalizes it, and applies RINT to the
expression values. The output is a pseudobulk expression matrix saved as a TSV file.    
"""

from scipy.stats import rankdata, norm
import scanpy as sc
import pandas as pd
import sys
import os
import numpy as np

# Perform Rank Inverse Normalization Transformation (RINT) with k=0.5
def RINT(x):
    ranks = rankdata(x, method='average')
    rint = norm.ppf((ranks - 0.5) / len(ranks))
    return rint

tissue = sys.argv[1]
ancestry = sys.argv[2]
celltype = sys.argv[3]

sample_filter = pd.read_csv(f"/path/to/{tissue}/sc-eQTL/genotype/{tissue}_{ancestry}_indivpruned.fam", sep='\t', header=None)[0].tolist()
adata = sc.read_h5ad(f'/path/to/{tissue}/h5ad_celltype/{tissue}_rawQC_{celltype}_{ancestry}.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata)

expr_out = f"/path/to/celltype_matrix/norm_rint_tsv/{tissue}_expr_norm_rint_{ancestry}_{celltype}.tsv"
if not os.path.isfile(expr_out):
    expression_df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    pseudobulk_df = expression_df.groupby(adata.obs['project_sample']).mean()
    pseudobulk_df = pseudobulk_df.loc[pseudobulk_df.index.isin(sample_filter),:] 
    pseudobulk_df_RINT = pseudobulk_df.apply(lambda x: RINT(x), axis=0)
    print("Write out normalized expression matrix")
    pseudobulk_df_RINT.T.to_csv(expr_out, sep='\t', index=True, header=True, index_label='gene')
