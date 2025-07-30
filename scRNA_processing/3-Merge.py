#!/usr/bin/env python3
"""
scRNA-seq Sample Merging Pipeline

This script combines multiple scRNA-seq sample AnnData objects into a single merged dataset.
It reads a list of sample paths from an input file and concatenates them with proper metadata.

Usage:
    python merge_samples.py <samples_list_file>

Example:
    python merge_samples.py normal_samples.txt
"""

import os
import sys
import time
from typing import List

import scanpy as sc
import anndata as ad
import pandas as pd
from progress.bar import Bar

def load_sample_paths(file_path: str) -> tuple:

    """
    Load sample file paths from a text file.
    
    Args:
        file_path: Path to file containing h5ad files to merge
        
    Returns:
        Tuple of (list of file paths, list of sample names)

    """
    h5ad_files = pd.read_csv(file_path, header=None)[0].tolist()
    sample_names = [os.path.basename(file).replace('_rawQC.h5ad', '') for file in h5ad_files]
    
    return h5ad_files, sample_names

def load_and_process_samples(sample_names: List[str], root_dir: str, progress_bar=None) -> List[ad.AnnData]:
    """
    Load and process AnnData objects for each sample.

    
    Args:
        sample_names: List of sample identifiers
        root_dir: Base directory containing h5ad files
        progress_bar: Optional progress tracking object
        
    Returns:
        List of processed AnnData objects
    """
    adatas = []
    
    for sample in sample_names:
        print(f"Processing {sample}...")
        path = f'{root_dir}/{sample}_rawQC.h5ad'
        
        # Load the sample data
        adata = sc.read_h5ad(path)
        
        # Keep only essential metadata columns
        adata.obs = adata.obs[[
            'project', 'sample', 'n_genes_by_counts', 
            'total_counts', 'pct_counts_mt', 'celltypist_pred'
        ]]
        
        # Add to collection
        adatas.append(adata)
        
        # Update progress bar
        if progress_bar:
            progress_bar.suffix = f'sample id: {sample}'
            progress_bar.next()
    
    return adatas

def main():
    """Main function to execute the merging workflow."""
    # Parse command line arguments
    if len(sys.argv) != 3
        print("Usage: python merge_samples.py <tissue> <samples_list_file>")
        sys.exit(1)
    
    tissue = sys.argv[1]
    merge_list = sys.argv[2]
    root_dir = f'/storage/yangjianLab/chenchang/scRNA/{tissue}/h5ad/'
    
    # Load sample information
    h5ad_files, sample_names = load_sample_paths(f'{root_dir}/{merge_list}')
    
    # Set up progress tracking
    progress_bar = Bar("Merge samples", max=len(sample_names))
    progress_bar.check_tty = False
    
    # Time the merging process
    start_time = time.time()
    
    # Load and process samples
    adatas = load_and_process_samples(sample_names, root_dir, progress_bar)

    
    # Perform the merge operation
    print(f"Merging {len(adatas)} samples...")
    merged_adata = ad.concat(adatas, join='outer', index_unique=None)
    
    # Save the merged result
    output_path = f'/storage/yangjianLab/chenchang/scRNA/{tissue}/h5ad_merge/{tissue}_merged_{merge_list}.h5ad'
    print(f"Writing merged data to {output_path}")
    merged_adata.write_h5ad(output_path)
    
    # Report timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time elapsed for merging: {elapsed_time:.2f} seconds")
    print(f"Merged data contains {merged_adata.n_obs} cells and {merged_adata.n_vars} genes")

if __name__ == "__main__":
    main()