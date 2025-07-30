#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Single-cell Variational Inference (scVI) and Single-cell Annotation
using Variational Inference (scANVI) implementation for integration,
analysis, and visualization of single-cell RNA-seq data.

This version allows specifying the tissue type as a command-line argument.
"""

import os
import sys
import shutil
import logging
import numpy as np
import scanpy as sc
import pandas as pd
import scvi

def load_and_prepare_data(tissue, hvg_suffix):
    """
    Load data from h5ad file and prepare it for analysis.


    Parameters
    ----------
    tissue : str
        The name of the tissue being analyzed (e.g., 'Colon').
    hvg_suffix : str
        Suffix for the highly variable genes file.


    Returns
    -------
    adata : AnnData
        Annotated data object with merged metadata.
    """
    print("Loading and preparing data...")

    file_path = f"/storage/yangjianLab/chenchang/scRNA/{tissue}/h5ad_merge/{tissue}_merged_rawQC_{hvg_suffix}.h5ad"
    adata = sc.read_h5ad(file_path)

    print(f'Number of cells: {adata.shape[0]}')

    return adata

def train_scvi_model(adata, batch_key, model_suffix):
    """
    Set up and train an scVI model on the data.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.
    batch_key : str
        Key in adata.obs for batch correction.
    model_suffix : str
        Suffix for model saving.

    Returns
    -------
    model : scvi.model.SCVI
        Trained scVI model.
    """
    if os.path.exists(f'../integration/model/scvi_{model_suffix}'):
        # Load the existing model
        print("Load existing scVI model")
        model = scvi.model.SCVI.load(f'../integration/model/scvi_{model_suffix}', adata)

    else:
        print("Setting up and training scVI model...")
        categorical_covariate_keys = ["platform"]
        continuous_covariate_keys = ['pct_counts_mt','total_counts']
        # Set up AnnData object for scVI
        scvi.model.SCVI.setup_anndata(
            adata,
            batch_key=batch_key,
            categorical_covariate_keys=categorical_covariate_keys,
            continuous_covariate_keys=continuous_covariate_keys
        )

        # Initialize model
        model = scvi.model.SCVI(
            adata,
            n_layers=2,
            n_latent=30,
            n_hidden=128,
            dispersion='gene',
            gene_likelihood="nb"
        )

        # Train model
        model.train(
            batch_size=512,
            max_epochs=500,
            plan_kwargs={"lr": 0.001, "lr_patience": 10, "reduce_lr_on_plateau": True},
            early_stopping=True,
            early_stopping_patience=10
        )

        print(model)
        print('Saving scVI model...')
        model.save(f'../integration/model/scvi_{model_suffix}')

    return model

def train_scanvi_model(model, adata, labels_key, model_suffix):
    """
    Set up and train an scANVI model from an existing scVI model.

    Parameters
    ----------
    model : scvi.model.SCVI
        Trained scVI model.
    adata : AnnData
        Annotated data object.
    labels_key : str
        Key in adata.obs containing cell type labels.
    model_suffix : str
        Suffix for model saving.

    Returns
    -------
    str
        The name of the saved model file.
    """
    print("Setting up and training scANVI model...")

    # Initialize scANVI model from scVI model
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model,
        adata=adata,
        labels_key=labels_key,
        unlabeled_category="Unknown"
    )

    # Set up checkpoint callback
    checkpoint_dir = "../integration/checkpoints"
    checkpoint_filename = f"scanvi_{model_suffix}_checkpoint"
    checkpoint_callback = scvi.train.SaveCheckpoint(
        dirpath=checkpoint_dir,
        filename=checkpoint_filename,
        monitor="validation_accuracy",
        load_best_on_end=True,
        mode="max"
    )

    model_name = f'scanvi_{model_suffix}'

    try:
        # Train scANVI model
        scanvi_model.train(
            max_epochs=200,
            plan_kwargs={"weight_decay": 0, "lr": 0.0001, "reduce_lr_on_plateau": True, "lr_patience": 10},
            early_stopping_monitor="validation_accuracy",
            early_stopping=True,
            early_stopping_patience=10,
            enable_checkpointing=True,
            early_stopping_mode='max',
            callbacks=[checkpoint_callback]
        )

        # Save model and embeddings
        scanvi_model.save(f'../integration/model/{model_name}')
        scanvi_embedding = scanvi_model.get_latent_representation(adata)
        scanvi_embedding = pd.DataFrame(scanvi_embedding, index=adata.obs_names)
        scanvi_embedding.to_csv(
            f"../integration/embedding/{model_name}.csv.gz",
            header=None,
            compression="gzip"
        )
        print(scanvi_model)

    except Exception as e:
        print(f"Error during scANVI training: {e}")
        # Extract the best model from checkpoints
        best_model_path = checkpoint_callback.best_model_path
        if os.path.exists(best_model_path):
            os.makedirs(os.path.dirname(f'../integration/model/{model_name}_best'), exist_ok=True)
            shutil.move(best_model_path, f'../integration/model/{model_name}_best')
            model_name = f"{model_name}_best"
            print(f"Best model moved to ../integration/model/{model_name}")
        # Load best model and generate embeddings
        best_scanvi_model = scvi.model.SCANVI.load(
            f'../integration/model/{model_name}',
            adata=adata
        )

        best_scanvi_embedding = best_scanvi_model.get_latent_representation(adata)
        best_scanvi_embedding = pd.DataFrame(best_scanvi_embedding, index=adata.obs_names)
        best_scanvi_embedding.to_csv(
            f"../integration/embedding/{model_name}.csv.gz",
            header=None,
            compression="gzip"
        )
        print(f"Embeddings for best model saved as {model_name}.csv.gz")

    return model_name

def create_umap_and_clustering(tissue, model_name, hvg_suffix, n_neighbors, leiden_res, adata=None):
    """
    Generate UMAP visualization and Leiden clustering from scANVI embeddings.


    Parameters
    ----------
    tissue : str
        The name of the tissue being analyzed (e.g., 'Colon').
    model_name : str
        Name of the model embeddings to use.
    hvg_suffix : str
        Suffix for the highly variable genes file.
    n_neighbors : int
        Number of neighbors for nearest neighbors calculation.
    leiden_res : float
        Resolution parameter for Leiden clustering.
    adata : AnnData, optional
        If provided, uses this AnnData object instead of loading from disk.
    """
    # Set up logging
    os.makedirs("job_reports", exist_ok=True)
    logging.basicConfig(filename=f"job_reports/UMAP_{model_name}.log",
                        filemode='a',
                        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.INFO)

    logging.info('Reading data')

    h5ad_nn_path = f'../h5ad_merge/{tissue}_merged_rawQC_{model_name}_nn{n_neighbors}.h5ad'
    h5ad_raw_path = f"../h5ad_merge/{tissue}_merged_rawQC_{hvg_suffix}.h5ad"

    # Check if we need to load data
    if adata is None:
        # Check if processed data already exists
        if os.path.isfile(h5ad_nn_path):
            adata = sc.read_h5ad(h5ad_nn_path)
        else:
            # Load embeddings and original data
            scanvi_embedding = pd.read_csv(f"../integration/embedding/{model_name}.csv.gz", header=None, index_col=0)
            adata = sc.read_h5ad(h5ad_raw_path)
            adata = adata[scanvi_embedding.index].copy()
            adata.obsm["X_scANVI"] = np.array(scanvi_embedding)
            # Calculate nearest neighbors
            logging.info('Calculating nearest neighbors')
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_scANVI', key_added="scanvi")

            # Save intermediate result
            adata.write_h5ad(h5ad_nn_path)
    else:
        # Use provided adata but make sure it has the scANVI embeddings and neighbors
        if "X_scANVI" not in adata.obsm:
            scanvi_embedding = pd.read_csv(f"../integration/embedding/{model_name}.csv.gz", header=None, index_col=0)
            adata = adata[scanvi_embedding.index].copy()
            adata.obsm["X_scANVI"] = np.array(scanvi_embedding)

        # Calculate nearest neighbors if not already done
        if "scanvi" not in adata.uns:
            logging.info('Calculating nearest neighbors')
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_scANVI', key_added="scanvi")

            # Save intermediate result
            adata.write_h5ad(h5ad_nn_path)


    # Generate UMAP if not already done
    logging.info('Generating UMAP')
    umap_path = f"../integration/umap/{model_name}_nn{n_neighbors}_umap.csv.gz"
    if os.path.isfile(umap_path) and "X_umap" not in adata.obsm:
        umap = pd.read_csv(umap_path, header=None, index_col=0)
        adata.obsm["X_umap"] = np.array(umap)
    elif "X_umap" not in adata.obsm:
        sc.tl.umap(adata, neighbors_key="scanvi", min_dist=0.3)
        os.makedirs("../integration/umap", exist_ok=True)
        pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names).to_csv(
            umap_path,
            compression="gzip",
            header=None
        )

    # Perform Leiden clustering
    logging.info('Performing Leiden clustering')
    leiden_key = f"leiden_{leiden_res}"
    if leiden_key not in adata.obs:
        sc.tl.leiden(
            adata,
            resolution=leiden_res,
            neighbors_key="scanvi",
            key_added=leiden_key,
            flavor='igraph',
            n_iterations=2
        )

    # Save clustering results
    os.makedirs("../integration/leiden", exist_ok=True)
    adata.obs[leiden_key].to_csv(
        f'../integration/leiden/leiden_{leiden_res}_{model_name}.csv.gz',
        compression="gzip"
    )

    # Save complete data
    logging.info('Saving final data')
    adata.write_h5ad(f'../h5ad_merge/{tissue}_merged_rawQC_{model_name}_umap_leiden.h5ad')
    # Generate plots
    logging.info("Generating plots")
    os.makedirs("../integration/figures", exist_ok=True)
    sc.settings.figdir = "../integration/figures/"
   
    sc.pl.umap(
        adata,
        color=[leiden_key, 'celltypist_level3'],
        save=f"_{tissue}_{model_name}_celltypist.png"
    )

def run_full_workflow(tissue, hvg_suffix, batch_key, model_suffix, n_neighbors, leiden_res):
    """
    Run the complete workflow from training scVI/scANVI to generating UMAP and clustering.

    Parameters
    ----------
    tissue : str
        The name of the tissue being analyzed (e.g., 'Colon').
    hvg_suffix : str
        Suffix for the highly variable genes file.
    batch_key : str
        Key in adata.obs for batch correction.
    model_suffix : str
        Suffix for model saving.
    n_neighbors : int
        Number of neighbors for nearest neighbors calculation.
    leiden_res : float
        Resolution parameter for Leiden clustering.
    """
    # Create necessary directories
    os.makedirs("../integration/model", exist_ok=True)
    os.makedirs("../integration/embedding", exist_ok=True)
    os.makedirs("../integration/checkpoints", exist_ok=True)

    # Load and prepare data
    adata = load_and_prepare_data(tissue, hvg_suffix)
    # Train scVI model
    model = train_scvi_model(adata, batch_key, model_suffix)

    # Train scANVI model
    model_name = train_scanvi_model(model, adata, labels_key="celltypist_level3", model_suffix=model_suffix)

    # Get scANVI embeddings and add to AnnData
    scanvi_embedding = pd.read_csv(f"../integration/embedding/{model_name}.csv.gz",
                                   header=None, index_col=0)
    adata = adata[scanvi_embedding.index].copy()
    adata.obsm["X_scANVI"] = np.array(scanvi_embedding)

    # Generate UMAP and clustering, passing the already loaded AnnData
    create_umap_and_clustering(tissue, model_name, hvg_suffix, n_neighbors, leiden_res, adata=adata)

    print("Full workflow completed successfully!")

def main():
    """Main function to execute the scVI/scANVI workflow."""
    if len(sys.argv) < 2:
        print("Usage: ")
        print("  For scVI/scANVI training: python script.py train <tissue> <hvg_suffix> <batch_key> <model_suffix>")
        print("  For UMAP visualization: python script.py umap <tissue> <model_name> <hvg_suffix> <n_neighbors> <leiden_res>")
        print("  For complete workflow: python script.py full <tissue> <hvg_suffix> <batch_key> <model_suffix> <n_neighbors> <leiden_res>")
        sys.exit(1)

    mode = sys.argv[1]


    if mode == "train":
        if len(sys.argv) != 6:
            print("Usage for training mode: python script.py train <tissue> <hvg_suffix> <batch_key> <model_suffix>")
            sys.exit(1)

        tissue = sys.argv[2]
        hvg_suffix = sys.argv[3]
        batch_key = sys.argv[4]  # project or project_sample
        model_suffix = sys.argv[5]


        # Create necessary directories
        os.makedirs("../integration/model", exist_ok=True)
        os.makedirs("../integration/embedding", exist_ok=True)
        os.makedirs("../integration/checkpoints", exist_ok=True)
        # Load and prepare data
        adata = load_and_prepare_data(tissue, hvg_suffix)
        # Train scVI model
        model = train_scvi_model(adata, batch_key, model_suffix)
        # Train scANVI model
        train_scanvi_model(model, adata, labels_key="celltypist_level3", model_suffix=model_suffix)
    elif mode == "umap":
        if len(sys.argv) != 7:
            print("Usage for UMAP mode: python script.py umap <tissue> <model_name> <hvg_suffix> <n_neighbors> <leiden_res>")
            sys.exit(1)

        tissue = sys.argv[2]
        model_name = sys.argv[3]
        hvg_suffix = sys.argv[4]
        n_neighbors = int(sys.argv[5])
        leiden_res = float(sys.argv[6])


        # Create UMAP visualization and clustering
        create_umap_and_clustering(tissue, model_name, hvg_suffix, n_neighbors, leiden_res)


    elif mode == "full":
        if len(sys.argv) != 8:
            print("Usage for full workflow mode: python script.py full <tissue> <hvg_suffix> <batch_key> <model_suffix> <n_neighbors> <leiden_res>")
            sys.exit(1)

        tissue = sys.argv[2]
        hvg_suffix = sys.argv[3]
        batch_key = sys.argv[4]
        model_suffix = sys.argv[5]
        n_neighbors = int(sys.argv[6])
        leiden_res = float(sys.argv[7])

        # Run complete workflow
        run_full_workflow(tissue, hvg_suffix, batch_key, model_suffix, n_neighbors, leiden_res)


    else:
        print(f"Unknown mode: {mode}")
        print("Available modes: train, umap, full")
        sys.exit(1)

if __name__ == "__main__":
    main()
