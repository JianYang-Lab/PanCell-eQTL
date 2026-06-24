#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Level 1: Coarse Single-cell Variational Inference (scVI) and scANVI.
Integrates data, predicts coarse lineages using KNN, and generates UMAPs.
"""

import os
import sys
import shutil
import numpy as np
import scanpy as sc
import pandas as pd
import scvi
from sklearn.neighbors import KNeighborsClassifier

def load_and_prepare_data(tissue, hvg_suffix):
    print("Loading and preparing data...")
    file_path = f"{tissue}/h5ad_merge/{tissue}_merged_rawQC_{hvg_suffix}.h5ad"
    adata = sc.read_h5ad(file_path)
    
    print(f'Number of cells: {adata.shape[0]}')
    return adata

def mask_low_confidence_labels(adata, labels_key, conf_key, threshold=0.8):
    print(f"Masking labels in {labels_key} with confidence < {threshold} based on {conf_key}...")
    scanvi_label_key = f"{labels_key}_scanvi_input"
    adata.obs[scanvi_label_key] = adata.obs[labels_key].astype(str) 
    mask = adata.obs[conf_key] < threshold
    adata.obs.loc[mask, scanvi_label_key] = 'Unknown'
    adata.obs[scanvi_label_key] = adata.obs[scanvi_label_key].astype('category')
    print(f"Cells marked as Unknown: {mask.sum()} / {len(adata)}")
    return scanvi_label_key

def train_scvi_model(adata, tissue, batch_key, model_suffix):
    model_path = f'{tissue}/integration/model/scvi_{model_suffix}'
    if os.path.exists(model_path):
        print("Load existing scVI model")
        model = scvi.model.SCVI.load(model_path, adata)
    else:
        print("Setting up and training scVI model...")
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
            early_stopping=True, early_stopping_patience=10
        )
        print('Saving scVI model...')
        model.save(model_path)
    return model

def train_scanvi_model(model, adata, tissue, labels_key, model_suffix):
    model_path = f'{tissue}/integration/model/scanvi_{model_suffix}'
    model_name = f'scanvi_{model_suffix}'
    if os.path.exists(model_path):
        print("Load existing scANVI model")
    #    adata = sc.read_h5ad(f"{tissue}/h5ad_merge/{tissue}_merged_rawQC_scanvi_{model_suffix}_umap_leiden.h5ad")
        scanvi_model = scvi.model.SCANVI.load(model_path, adata)
    else:
        print("Setting up and training scANVI model...")
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            model, adata=adata, labels_key=labels_key, unlabeled_category="Unknown"
        )

        checkpoint_dir = f"{tissue}/integration/checkpoints"
        checkpoint_filename = f"scanvi_{model_suffix}_checkpoint"
        checkpoint_callback = scvi.train.SaveCheckpoint(
            dirpath=checkpoint_dir, filename=checkpoint_filename,
            monitor="validation_accuracy", load_best_on_end=True, mode="max"
        )

        try:
            scanvi_model.train(
                max_epochs=200,
                plan_kwargs={"weight_decay": 0, "lr": 0.0001, "reduce_lr_on_plateau": True, "lr_patience": 10},
                early_stopping_monitor="validation_accuracy", early_stopping=True,
                early_stopping_patience=10, enable_checkpointing=True,
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
        scanvi_embedding = pd.DataFrame(scanvi_embedding, index=adata.obs_names)
        scanvi_embedding.to_csv(f"{tissue}/integration/embedding/{model_name}.csv.gz", header=None, compression="gzip")
        
    return model_name

def run_knn_rescue(adata, tissue, latent_key, input_label_key, final_label_key, model_suffix):
    print("Running KNN to rescue Unknown cells...")
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
    adata.obs[final_label_key].to_csv(f"{tissue}/h5ad_merge/{tissue}_scanvi_{model_suffix}_celltypist_coarse_knn_final.txt",sep="\t")
    print("KNN rescue complete.")

def generate_umap_and_plot(adata, tissue, n_neighbors, leiden_res, color_keys, save_prefix):
    """Calculates neighbors, UMAP, Leiden clustering, and saves plots."""
    print(f"Calculating nearest neighbors (n={n_neighbors}) on scANVI latent space...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_scANVI', key_added="scanvi")

    print("Generating UMAP...")
    sc.tl.umap(adata, neighbors_key="scanvi", min_dist=0.3)

    print(f"Performing Leiden clustering (resolution={leiden_res})...")
    leiden_key = f"leiden_{leiden_res}"
    sc.tl.leiden(adata, resolution=leiden_res, neighbors_key="scanvi", key_added=leiden_key, flavor='igraph', n_iterations=2)
    
    # Add leiden to the plotting list
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

def run_full_workflow(tissue, hvg_suffix, batch_key, labels_key, conf_key, threshold, model_suffix, n_neighbors, leiden_res):
    os.makedirs(f"{tissue}/integration/model", exist_ok=True)
    os.makedirs(f"{tissue}/integration/embedding", exist_ok=True)
    os.makedirs(f"{tissue}/integration/checkpoints", exist_ok=True)

    adata = load_and_prepare_data(tissue, hvg_suffix)
    scanvi_input_key = mask_low_confidence_labels(adata, labels_key, conf_key, threshold)
    
    base_model_name = f"scanvi_{model_suffix}"
    best_model_name = f"scanvi_{model_suffix}_best"
    embedding_path_base = f"{tissue}/integration/embedding/{base_model_name}.csv.gz"
    embedding_path_best = f"{tissue}/integration/embedding/{best_model_name}.csv.gz"

    if os.path.isfile(embedding_path_best):
        print(f"Found existing best embedding: {embedding_path_best}")
        model_name = best_model_name
        embedding_path = embedding_path_best
    elif os.path.isfile(embedding_path_base):
        print(f"Found existing base embedding: {embedding_path_base}")
        model_name = base_model_name
        embedding_path = embedding_path_base
    else:
        print("No existing embedding found. Initiating training...")
        model = train_scvi_model(adata, tissue, batch_key, model_suffix)
        model_name = train_scanvi_model(model, adata, tissue, labels_key=scanvi_input_key, model_suffix=model_suffix)
        embedding_path = f"{tissue}/integration/embedding/{model_name}.csv.gz"
        
    scanvi_embedding = pd.read_csv(embedding_path, header=None, index_col=0)
    adata = adata[scanvi_embedding.index].copy()
    adata.obsm["X_scANVI"] = np.array(scanvi_embedding)

    final_coarse_key = f"{labels_key}_knn_final"
    run_knn_rescue(adata, tissue, "X_scANVI", scanvi_input_key, final_coarse_key, model_suffix)

    save_prefix = f"{tissue}_{model_name}_coarse"
    color_keys = [batch_key, scanvi_input_key, final_coarse_key]
    generate_umap_and_plot(adata, tissue, n_neighbors, leiden_res, color_keys, save_prefix)

    output_path = f'{tissue}/h5ad_merge/{tissue}_merged_{model_name}_coarse_annotated.h5ad'
    print(f"Saving coarse annotated data to {output_path}")
    adata.write_h5ad(output_path)
    print("Level 1 workflow completed successfully!")

def main():
    os.chdir("/path/to/your/directory")
    if len(sys.argv) != 10:
        print("Usage: python SCANVI_coarse.py <tissue> <hvg_suffix> <batch_key> <coarse_labels_key> <conf_score_key> <threshold> <model_suffix> <n_neighbors> <leiden_res>")
        sys.exit(1)

    tissue = sys.argv[1]
    hvg_suffix = sys.argv[2]
    batch_key = sys.argv[3]
    labels_key = sys.argv[4]
    conf_key = sys.argv[5]
    threshold = float(sys.argv[6])
    model_suffix = sys.argv[7]
    n_neighbors = int(sys.argv[8])
    leiden_res = float(sys.argv[9])

    run_full_workflow(tissue, hvg_suffix, batch_key, labels_key, conf_key, threshold, model_suffix, n_neighbors, leiden_res)

if __name__ == "__main__":
    main()