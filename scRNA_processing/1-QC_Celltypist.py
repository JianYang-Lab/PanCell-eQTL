#!/usr/bin/env python3

"""
scRNA-seq Analysis Pipeline: Quality Control and Cell Type Annotation

This script performs two main steps:
1) Quality control for scRNA count matrix per sample
2) Automatic cell type annotation using Celltypist

Usage:
    python script.py <tissue> <platform> [additional arguments based on platform]
"""

import os
import sys
import numpy as np

import pandas as pd
import scanpy as sc

import scrublet as scr
from scipy.sparse import csr_matrix

from scipy.stats import median_abs_deviation
import celltypist

from celltypist import models

def is_outlier(adata, metric, nmads):
    """Detect outliers based on median absolute deviation."""

    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (

        np.median(M) + nmads * median_abs_deviation(M) < M

    )
    return outlier

def calculate_expected_doublet_rate(ncell):
    """Calculate expected doublet rate based on number of cells."""

    if ncell < 500:

        return 0.004

    elif ncell < 1000:
        return 0.004

    elif ncell < 2000:
        return 0.008
    elif ncell < 3000:
        return 0.016

    elif ncell < 4000:
        return 0.023

    elif ncell < 5000:
        return 0.031
    elif ncell < 6000:

        return 0.039
    elif ncell < 7000:

        return 0.046

    elif ncell < 8000:
        return 0.054
    elif ncell < 9000:

        return 0.061
    elif ncell < 10000:

        return 0.069

    else:
        return 0.076

def load_data(platform, *args):
    """Load data based on the specified platform."""
    if platform == "10X":
        project, sample = args
        path_filter = os.path.join(project, sample, 'outs/filtered_feature_bc_matrix')

        adata = sc.read_10x_mtx(path_filter)
        
    elif platform == "10X_h5":

        project, sample = args

        path_filter = os.path.join(project, f'raw_matrix/{project}_{sample}.h5')

        adata = sc.read_10x_h5(path_filter)

        adata.var_names_make_unique()

        
    elif platform == "h5ad":
        project, sample = args

        path_filter = os.path.join(project, f'raw_matrix/{project}_{sample}.h5ad')
        adata = sc.read_h5ad(path_filter)
        adata.raw = None
        if 'orig.ident' in adata.obs.columns:

            adata.obs.drop(['orig.ident', 'nCount_RNA', 'nFeature_RNA'], axis=1, inplace=True)

        if 'features' in adata.var.columns:

            adata.var.drop('features', axis=1, inplace=True)

        
    elif platform == "10X_multi":

        project, multi, sample = args
        path_filter = os.path.join(
            project, multi, "outs/per_sample_outs", sample, 

            "count/sample_filtered_feature_bc_matrix"
        )
        adata = sc.read_10x_mtx(path_filter, gex_only=True)

        sample = f"{multi}_{sample}"
        
    elif platform == "10X_csv":

        import anndata as ad
        project, sample = args

        path_filter = f'{project}/raw_matrix/{project}_{sample}.counts.csv.gz'
        data_df = pd.read_csv(path_filter, index_col=0)

        sparse_matrix = csr_matrix(data_df)

        adata = ad.AnnData(X=sparse_matrix.transpose())
        adata.var_names = data_df.index

        adata.obs_names = data_df.columns
        
    elif platform == "aggr":
        project, sample = args

        path_filter = os.path.join(project, sample, "outs/count/filtered_feature_bc_matrix")
        adata = sc.read_10x_mtx(path_filter)

        
    elif platform == "multiaggr":
        project, sample = args
        path_filter = os.path.join(project, sample, "outs/count/filtered_feature_bc_matrix")
        adata = sc.read_10x_mtx(path_filter, gex_only=True)

        
    elif platform == "freemuxlet":
        project, multi, sample = args

        path_filter = os.path.join(project, multi, "outs/filtered_feature_bc_matrix")

        adata = sc.read_10x_mtx(path_filter)
        path_freemuxlet = os.path.join(project, multi, "freemuxlet", f"barcode_sample_{sample}.csv")
        barcodes = pd.read_csv(path_freemuxlet, header=None)[0]

        adata = adata[barcodes, :]
        sample = f"{multi}_{sample}"
        
    elif platform == "singleron":
        project, sample = args

        path_filter = os.path.join(project, sample, 'outs/filtered')

        adata = sc.read_10x_mtx(path_filter)

        
    elif platform == "smartseq":
        import anndata as ad

        project, sample = args
        path_filter = os.path.join(project, "raw_matrix", f"{project}_{sample}.genes.results")

        data_df = pd.read_csv(path_filter, sep="\t", index_col=0)
        sparse_matrix = csr_matrix(data_df)
        adata = ad.AnnData(X=sparse_matrix.transpose())

        adata.var_names = data_df.index
        adata.obs_names = data_df.columns    
        adata.var_names_make_unique()
        
    elif platform in ["BD", "dropseq"]:
        import anndata as ad
        project, sample = args

        separator = "\t"
        file_suffix = "counts.tsv.gz" if platform == "BD" else "dge.txt.gz"
        path_filter = f'{project}/raw_matrix/{project}_{sample}.{file_suffix}'
        data_df = pd.read_csv(path_filter, sep=separator, index_col=0)

        sparse_matrix = csr_matrix(data_df)
        adata = ad.AnnData(X=sparse_matrix.transpose())
        adata.var_names = data_df.index

        adata.obs_names = data_df.columns
        
    else:
        raise ValueError(f"Platform '{platform}' not recognized.")

    
    # Update cell barcodes with project and sample information
    adata.obs_names = [f"{name}-{project}-{sample}" for name in adata.obs_names]

    adata.obs["project"] = project

    adata.obs["sample"] = sample

    
    return adata, project, sample

def main():
    """Main function to run the scRNA-seq analysis pipeline."""

    # Parse command line arguments
    tissue = sys.argv[1]

    platform = sys.argv[2]
    
    # Set working directory
    os.chdir(f"/storage/yangjianLab/chenchang/scRNA/{tissue}/raw_data")

    
    print("Loading data...")

    adata, project, sample = load_data(platform, *sys.argv[3:])

    print(f"Platform: {platform}; Project: {project}; Sample: {sample}")

    
    # Basic filtering

    print("Performing quality control filtering...")
    sc.pp.filter_cells(adata, min_genes=200)

    sc.pp.filter_genes(adata, min_cells=3)

    
    # Annotate mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    #adata.var["rb"] = adata.var_names.str.startswith(("RPS","RPL"))
    #adata.var["hb"] = adata.var_names.str.startswith(("^HB[^(P)]"))
    
    sc.pp.calculate_qc_metrics(

        adata, qc_vars=["mt"], percent_top=[20], log1p=True, inplace=True

    )
    
    # Detect outliers

    adata.obs["outlier"] = (

        is_outlier(adata, "log1p_total_counts", 5) | 
        is_outlier(adata, "log1p_n_genes_by_counts", 5)
    )
    adata.obs.outlier.value_counts()

    
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 5) | (

        adata.obs["pct_counts_mt"] > 20
    )
    adata.obs.mt_outlier.value_counts()
    
    print(f"Total number of cells: {adata.n_obs}")

    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    print(f"Number of cells after filtering low-quality cells: {adata.n_obs}")
    
    # Doublet detection using scrublet
    print("Detecting doublets...")
    ncell = adata.shape[0]
    dbrate = calculate_expected_doublet_rate(ncell)
    
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=dbrate)

    doublet_scores, predicted_doublets = scrub.scrub_doublets(

        min_counts=2, 
        min_cells=3, 
        min_gene_variability_pctl=85, 
        n_prin_comps=30
    )
    
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["doublet_class"] = predicted_doublets

    adata = adata[~adata.obs.doublet_class].copy()
    
    ncell = adata.shape[0]
    if ncell >= 50:
        print(f'After filtering, {ncell} cells remaining.')

    else:
        raise ValueError("Number of cells is less than 50, discard the sample.")
    
    # Automated cell type annotation using CellTypist
    print("Performing automated cell type annotation using CellTypist...")
    adata_celltypist = adata.copy()

    sc.pp.normalize_total(adata_celltypist, target_sum=1e4)
    sc.pp.log1p(adata_celltypist)

    
    # Make .X dense for compatibility with celltypist
    adata_celltypist.X = adata_celltypist.X.toarray()

    
    # Load appropriate model for tissue
    #models.download_models(model=["Human_Lung_Atlas"])
    model = models.Model.load(model="Human_Lung_Atlas.pkl")
    predictions = celltypist.annotate(adata_celltypist, model=model, majority_voting=True)
    predictions_adata = predictions.to_adata()
    adata.obs["celltypist_pred"] = predictions_adata.obs.loc[adata.obs.index, "majority_voting"]
    
    # Save results
    print("Writing output file...")
    sample_h5ad_out = f"/storage/yangjianLab/chenchang/scRNA/{tissue}/h5ad/{project}_{sample}_rawQC.h5ad"

    adata.write_h5ad(sample_h5ad_out)
    
    print("Session finished")

if __name__ == "__main__":
    main()