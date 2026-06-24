#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Level 2: Fine Single-cell Variational Inference (scVI) and scANVI.
Operates on a specific subset, resolves high-resolution states, and generates subset UMAPs.
"""

import os
import sys
import gc
import shutil
import numpy as np
import scanpy as sc
import pandas as pd
import scvi
from sklearn.neighbors import KNeighborsClassifier
import celltypist
from celltypist import models
from joblib import Parallel, delayed

def load_and_subset_data(tissue, coarse_col, target_group, model_suffix):
    print(f"Loading and subsetting data for target group: {target_group}...")
    file_hvg_path = f"{tissue}/h5ad_merge/{tissue}_scanvi_{model_suffix}_merged_{target_group}_hvg.h5ad"
    
    if os.path.isfile(file_hvg_path):
        adata_sub = sc.read_h5ad(file_hvg_path)
    else:
        file_path = f"{tissue}/h5ad_merge/{tissue}_merged_new.h5ad"
        adata = sc.read_h5ad(file_path)
        
        if coarse_col not in adata.obs.columns:
            ct_coarse = pd.read_csv(f"{tissue}/h5ad_merge/{tissue}_scanvi_{model_suffix}_celltypist_coarse_knn_final.txt",sep="\t",index_col=0)
            adata.obs = pd.merge(adata.obs, ct_coarse, how="left", left_index=True, right_index=True) 
        
        if "platform" not in adata.obs.columns:
            platform = pd.read_csv(f"{tissue}/h5ad_merge/{tissue}_barcode_platform.txt", sep="\t", index_col=0)
            adata.obs = pd.merge(adata.obs, platform, how="left", left_index=True, right_index=True)
        
        adata_sub = adata[adata.obs[coarse_col] == target_group].copy()
        
        print("Fetch CellTypist fine labels...")
        ct_fine = celltypist_fine(tissue, adata_sub, target_group, model_suffix)
        adata_sub.obs = pd.merge(adata_sub.obs, ct_fine, how='left', left_index=True, right_index=True)
        
        print("Select HVGs...") 
        adata_sub = select_hvg(adata_sub, tissue, target_group, model_suffix)
        
    print(f'Number of cells in {target_group} subset: {adata_sub.shape[0]}')
    return adata_sub

def annotate_sample(adata_subset, model, prefix):
    """Preprocesses a single sample and runs all required celltypist models."""
    # Work on a copy to prevent modifying the original slice
    adata_ct = adata_subset.copy()
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct) 
    
    if hasattr(adata_ct.X, "toarray"):
        adata_ct.X = adata_ct.X.toarray()
    
    predictions = celltypist.annotate(adata_ct, model=model, majority_voting=False) 
    predictions_adata = predictions.to_adata(prefix=prefix)

    col_keep = [f"{prefix}predicted_labels", f"{prefix}conf_score"]

    return predictions_adata.obs[col_keep] 

def celltypist_fine(tissue, adata_sub, target_group, model_suffix):
    ct_fine_path = f"{tissue}/h5ad_merge/{tissue}_scanvi_{model_suffix}_celltypist_fine_{target_group}.txt"
    if not os.path.isfile(ct_fine_path):
        print(f"Predict fine labels within {target_group}...")
        
        model = models.Model.load(f"celltypist_model/{tissue}/{tissue}_{target_group}.pkl")
        
        prefix = "celltypist_fine_"
        samples_to_process = adata_sub.obs["project_sample"].unique() # Remove [0:2] to run all
        print(f"Running Celltypist on {len(samples_to_process)} samples in parallel...")

        results_list = Parallel(n_jobs=4)(
            delayed(annotate_sample)(adata_sub[adata_sub.obs["project_sample"] == ps], model, prefix) 
            for ps in samples_to_process
        )
        
        ct_fine = pd.concat(results_list, axis=0) 
        ct_fine.to_csv(ct_fine_path, sep="\t")
    else:
        print("CellTypist fine label already exists...")
        ct_fine = pd.read_csv(ct_fine_path, sep="\t", index_col=0)
    
    return ct_fine
    
def select_hvg(adata_sub, tissue, target_group, model_suffix):
    batch_key="project"
    file_hvg_path = f"{tissue}/h5ad_merge/{tissue}_scanvi_{model_suffix}_merged_{target_group}_hvg.h5ad"
    
    mt_genes = list(adata_sub.var_names[adata_sub.var_names.str.match(r'^MT-')])
    rp_genes = list(adata_sub.var_names[adata_sub.var_names.str.match(r'^RP[SL]')]) 
    ncRNA_genes = list(adata_sub.var_names[adata_sub.var_names.str.match(r'^[A-Z][A-Z][0-9].*\.[0-9]')])
    LINC_genes = list(adata_sub.var_names[adata_sub.var_names.str.match(r'(^LOC|LINC)[1-9]*')])
    remove_genes = mt_genes + rp_genes + ncRNA_genes + LINC_genes
   
    # Remove unused categories
    if adata_sub.obs[batch_key].dtype.name == 'category':
        adata_sub.obs[batch_key] = adata_sub.obs[batch_key].cat.remove_unused_categories()
    else:
        adata_sub.obs[batch_key] = adata_sub.obs[batch_key].astype('category')
        
    adata_ct = adata_sub.copy()
    
    # Filter out batches with too few cells (Scanpy HVG needs >= 2, scVI prefers more) 
    batch_counts = adata_ct.obs[batch_key].value_counts()
    valid_batches = batch_counts[batch_counts >= 3].index
    if len(valid_batches) < len(batch_counts):
        dropped_batches = list(set(batch_counts.index) - set(valid_batches))
        print(f"Warning: Dropping batches with <3 cells to prevent HVG crash: {dropped_batches}")
        adata_ct = adata_ct[adata_ct.obs[batch_key].isin(valid_batches)].copy()
        adata_ct.obs[batch_key] = adata_ct.obs[batch_key].cat.remove_unused_categories()
    
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct) 
    
    sc.pp.highly_variable_genes(adata_ct, flavor='seurat', batch_key="project", n_top_genes=2500)
    hvg_list = adata_ct.var_names[adata_ct.var["highly_variable"]] 
    hvg_keep = [g for g in hvg_list if g not in remove_genes]
    
    pd.Series(hvg_list).to_csv(f"{tissue}/h5ad_merge/{tissue}_scanvi_{model_suffix}_{target_group}_hvg_all.txt", sep="\t", index=None)
    pd.Series(hvg_keep).to_csv(f"{tissue}/h5ad_merge/{tissue}_scanvi_{model_suffix}_{target_group}_hvg_keep.txt", sep="\t", index=None)
    
    del adata_ct
    gc.collect() 
    
    adata_sub = adata_sub[:, adata_sub.var_names.isin(hvg_keep)].copy()
    adata_sub.write_h5ad(file_hvg_path)
    return adata_sub 

def mask_low_confidence_labels(adata, labels_key, conf_key, threshold=0.8):
    print(f"Masking labels in {labels_key} with confidence < {threshold}...")
    scanvi_label_key = f"{labels_key}_scanvi_input"
    adata.obs[scanvi_label_key] = adata.obs[labels_key].astype(str)
    
    mask = adata.obs[conf_key] < threshold
    adata.obs.loc[mask, scanvi_label_key] = 'Unknown'
    adata.obs[scanvi_label_key] = adata.obs[scanvi_label_key].astype('category')
    return scanvi_label_key

def train_scvi_model(adata, tissue, batch_key, model_suffix, target_group):
    model_path = f'{tissue}/integration/model/scvi_{model_suffix}_{target_group}'
    if os.path.exists(model_path):
        print(f"Load existing scVI model for {target_group}")
        model = scvi.model.SCVI.load(model_path, adata)
    else:
        print(f"Setting up and training scVI model for {target_group}...")
        categorical_covariate_keys = ["platform", "project"]
        continuous_covariate_keys = ['pct_counts_mt', "total_counts","n_genes_by_counts"] 
       
        categorical_covariate_keys = [k for k in categorical_covariate_keys if k in adata.obs.columns]
        continuous_covariate_keys = [k for k in continuous_covariate_keys if k in adata.obs.columns]
         
        scvi.model.SCVI.setup_anndata(
            adata, batch_key=batch_key,
            categorical_covariate_keys=categorical_covariate_keys,
            continuous_covariate_keys=continuous_covariate_keys
        )

        model = scvi.model.SCVI(
            adata, n_layers=2, n_latent=30, n_hidden=128,
            dispersion='gene', gene_likelihood="nb"
        )

        model.train(
            batch_size=128, max_epochs=500,
            plan_kwargs={"lr": 0.001, "lr_patience": 10, "reduce_lr_on_plateau": True},
            early_stopping=True, early_stopping_patience=15
        )
        model.save(model_path, overwrite=True)
    return model

def train_scanvi_model(model, adata, tissue, labels_key, model_suffix, target_group):
    model_path = f'{tissue}/integration/model/scanvi_{model_suffix}_{target_group}'
    model_name = f'scanvi_{model_suffix}_{target_group}'
    
    if os.path.exists(model_path):
        print(f"Loading existing scANVI model for {target_group}...")
        scanvi_model = scvi.model.SCANVI.load(model_path, adata)
    else:
        print(f"Setting up and training scANVI model for {target_group}...")
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            model, adata=adata, labels_key=labels_key, unlabeled_category="Unknown"
        )

        checkpoint_dir = f"{tissue}/integration/checkpoints"
        checkpoint_filename = f"scanvi_{model_suffix}_{target_group}_checkpoint"
        checkpoint_callback = scvi.train.SaveCheckpoint(
            dirpath=checkpoint_dir, filename=checkpoint_filename,
            monitor="validation_accuracy", load_best_on_end=True, mode="max"
        )

        try:
            scanvi_model.train(
                max_epochs=200,
                plan_kwargs={"weight_decay": 0, "lr": 0.0001, "reduce_lr_on_plateau": True, "lr_patience": 10},
                early_stopping_monitor="validation_accuracy", early_stopping=True,
                early_stopping_patience=15, enable_checkpointing=True,
                early_stopping_mode='max', callbacks=[checkpoint_callback]
            )
        except Exception as e:
            print(f"Error or early stop during scANVI: {e}")
            best_model_path = checkpoint_callback.best_model_path
            if os.path.exists(best_model_path):
                os.makedirs(os.path.dirname(f'{tissue}/integration/model/{model_name}_best'), exist_ok=True)
                shutil.move(best_model_path, f'{tissue}/integration/model/{model_name}_best')
                model_name = f"{model_name}_best"
            scanvi_model = scvi.model.SCANVI.load(f'{tissue}/integration/model/{model_name}', adata=adata)

        scanvi_model.save(f'{tissue}/integration/model/{model_name}', overwrite=True)
        
        scanvi_embedding = scanvi_model.get_latent_representation(adata)
        pd.DataFrame(scanvi_embedding, index=adata.obs_names).to_csv(
            f"{tissue}/integration/embedding/{model_name}.csv.gz", header=None, compression="gzip"
        )
        
    return scanvi_model, model_name

def run_knn_rescue(adata, tissue, latent_key, input_label_key, final_label_key, target_group, model_suffix):
    print("Running KNN to rescue Unknown fine cells...")
    known_mask = adata.obs[input_label_key] != 'Unknown'
    unknown_mask = adata.obs[input_label_key] == 'Unknown'
    
    X_known = adata.obsm[latent_key][known_mask]
    y_known = adata.obs[input_label_key][known_mask].values
    X_unknown = adata.obsm[latent_key][unknown_mask]
    
    adata.obs[final_label_key] = adata.obs[input_label_key].astype(str)
    
    if len(X_unknown) > 0 and len(X_known) > 0:
        n_neighbors = min(15, len(X_known))
        knn = KNeighborsClassifier(n_neighbors=n_neighbors, weights='distance', n_jobs=8)
        knn.fit(X_known, y_known)
        predicted_unknowns = knn.predict(X_unknown)
        adata.obs.loc[unknown_mask, final_label_key] = predicted_unknowns
        
    adata.obs[final_label_key] = adata.obs[final_label_key].astype('category')
    adata.obs[final_label_key].to_csv(f"{tissue}/h5ad_merge/{tissue}_scanvi_{model_suffix}_celltypist_fine_knn_final_{target_group}.txt",sep="\t")
    print("KNN fine rescue complete.")

def generate_umap_and_plot(adata, tissue, n_neighbors, leiden_res, color_keys, save_prefix):
    print(f"Calculating nearest neighbors (n={n_neighbors}) on scANVI latent space...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_scANVI', key_added="scanvi")

    print("Generating UMAP...")
    sc.tl.umap(adata, neighbors_key="scanvi", min_dist=0.3)

    print(f"Performing Leiden clustering (resolution={leiden_res})...")
    leiden_key = f"leiden_{leiden_res}"
    sc.tl.leiden(adata, resolution=leiden_res, neighbors_key="scanvi", key_added=leiden_key, flavor='igraph', n_iterations=2)
    
    if leiden_key not in color_keys:
        color_keys.append(leiden_key)

    print("Saving embeddings and plots...")
    os.makedirs(f"{tissue}/integration/umap", exist_ok=True)
    os.makedirs(f"{tissue}/integration/leiden", exist_ok=True)
    os.makedirs(f"{tissue}/integration/figures", exist_ok=True)

    pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names).to_csv(f"{tissue}/integration/umap/{save_prefix}_umap.csv.gz", compression="gzip", header=None)
    adata.obs[leiden_key].to_csv(f'{tissue}/integration/leiden/{save_prefix}_leiden.csv.gz', compression="gzip")

    sc.settings.figdir = f"{tissue}/integration/figures/"
    sc.pl.umap(adata, color=color_keys, save=f"_{save_prefix}.png", show=False)

def run_full_sub_workflow(tissue, coarse_col, target_group, batch_key, fine_labels_key, conf_key, threshold, model_suffix, n_neighbors, leiden_res):
    adata_sub = load_and_subset_data(tissue, coarse_col, target_group, model_suffix)
    
    if len(adata_sub) < 100:
        print(f"Subset {target_group} too small ({len(adata_sub)} cells). Skipping scVI integration.")
        return

    scanvi_input_key = mask_low_confidence_labels(adata_sub, fine_labels_key, conf_key, threshold)
    
    # Start Model Processing
    model = train_scvi_model(adata_sub, tissue, batch_key, model_suffix, target_group)
    scanvi_model, model_name = train_scanvi_model(model, adata_sub, tissue, labels_key=scanvi_input_key, model_suffix=model_suffix, target_group=target_group)

    print("Extracting scANVI latent representation...")
    scanvi_embedding = scanvi_model.get_latent_representation(adata_sub)
    adata_sub.obsm["X_scANVI"] = scanvi_embedding

    final_fine_key = f"{fine_labels_key}_knn_final"
    run_knn_rescue(adata_sub, tissue, "X_scANVI", scanvi_input_key, final_fine_key, target_group, model_suffix)

    save_prefix = f"{tissue}_{model_name}_{target_group}_fine"
    color_keys = [batch_key, scanvi_input_key, final_fine_key]
    generate_umap_and_plot(adata_sub, tissue, n_neighbors, leiden_res, color_keys, save_prefix)

    output_path = f'{tissue}/h5ad_merge/{tissue}_merged_{model_name}_fine_annotated.h5ad'
    print(f"Saving fine annotated subset data to {output_path}")
    adata_sub.write_h5ad(output_path)
    print(f"Level 2 workflow for {target_group} completed successfully!")

def main():
    os.chdir(f"/path/to/your/directory")
    if len(sys.argv) != 11:
        print("Usage: python SCANVI_fine.py <tissue> <coarse_col> <target_group> <batch_key> <fine_labels_key> <fine_conf_key> <threshold> <model_suffix> <n_neighbors> <leiden_res>")
        sys.exit(1)

    tissue = sys.argv[1]  
    coarse_col = sys.argv[2]           
    target_group = sys.argv[3]         
    batch_key = sys.argv[4]
    fine_labels_key = sys.argv[5]      
    fine_conf_key = sys.argv[6]
    threshold = float(sys.argv[7])
    model_suffix = sys.argv[8]         
    n_neighbors = int(sys.argv[9])
    leiden_res = float(sys.argv[10])

    run_full_sub_workflow(tissue, coarse_col, target_group, batch_key, fine_labels_key, fine_conf_key, threshold, model_suffix, n_neighbors, leiden_res)

if __name__ == "__main__":
    main()