import scanpy as sc 
import pandas as pd 
import numpy as np 

"""
Helper functions to assess annotation consistency and identify marker genes.
"""

def cluster_entropy(adata, leiden_label, leiden_value, label):
    """
    Calculates the normalized Shannon entropy for a specific cluster.
    """
    # 1. Determine the maximum possible entropy for normalization
    n_classes = adata.obs[label].nunique()
    
    # Handle edge case: if there's only 1 or 0 possible labels, entropy is 0
    if n_classes <= 1:
        return 0.0
        
    cell_count = sum(adata.obs[leiden_label] == leiden_value)
    
    # Handle edge case: empty cluster
    if cell_count == 0:
        return 0.0
        
    # 2. Calculate probabilities of each label within the cluster
    p = adata.obs.loc[adata.obs[leiden_label] == leiden_value, label].value_counts() / cell_count
    p = p[p > 0]
    
    # 3. Calculate raw Shannon entropy (using base 2)
    raw_entropy = -np.sum(p * np.log2(p))
    
    # 4. Normalize (Raw Entropy / Max Possible Entropy)
    max_entropy = np.log2(n_classes)
    normalized_entropy = raw_entropy / max_entropy
    
    return normalized_entropy


def get_high_entropy_clusters(adata, leiden_label, label, label_entropy_thres):
    # Calculate entropy for each new subcluster
    print(f"Calculate cluster entropy for {leiden_label}")
    label_entropy = {}
    
    for leiden_value in adata.obs[leiden_label].unique():
        label_entropy[leiden_value] = cluster_entropy(adata, leiden_label, leiden_value, label)
            
    high_entropy_clusters = [leiden for leiden, entropy in label_entropy.items() if entropy > label_entropy_thres]  
    return high_entropy_clusters


def recursive_subclustering(adata, leiden_values, level, n_neighbors, resolution, label_entropy_thres, parent_key, label_key):
    """
    Recursively perform subclustering on high-entropy clusters using a normalized threshold.
    
    Parameters:
    adata: AnnData object
    leiden_values: list of cluster IDs to subcluster
    level: current subclustering level (1-5)
    n_neighbors: number of neighbors for graph construction
    label_entropy_thres: normalized entropy threshold (0.0 to 1.0) for further subclustering
    parent_key: key for parent leiden clustering in adata.obs
    label_key: key for the predicted fine labels in adata.obs
    """
    if level > 5 or not leiden_values:
        return
    
    # Current key for storing subclusters
    current_key = f'leiden_l{level}'
    ld_sub_list_current = []

    print(f"\n--- Level {level} subclustering ---")

    # Perform subclustering on each cluster
    for leiden_value in leiden_values:
        print(f"Subclustering cluster: {leiden_value}")

        # Subset data for this cluster
        adata_sub = adata[adata.obs[parent_key]==leiden_value].copy()

        print("Computing neighbors")
        sc.pp.neighbors(adata_sub, n_neighbors=n_neighbors, use_rep='X_scANVI', key_added="scanvi")

        print("Running Leiden")
        sc.tl.leiden(adata_sub, 
                    neighbors_key="scanvi", 
                    resolution=resolution,
                    key_added="leiden_sub", 
                    flavor='igraph', 
                    n_iterations=2)

        # Create new subcluster labels
        ld_sub = pd.DataFrame({current_key: [f"{ld1}_{ld2}" for ld1, ld2 in 
                                        zip(adata_sub.obs[parent_key], adata_sub.obs['leiden_sub'])]})
        ld_sub.index = adata_sub.obs_names

        ld_sub_list_current.append(ld_sub)
    
    # Merge new labels into adata.obs
    if ld_sub_list_current:
        ld_sub = pd.concat(ld_sub_list_current)

        adata.obs = pd.merge(adata.obs, ld_sub, left_index=True, right_index=True, how='left')
        adata.obs.loc[adata.obs[current_key].isna(), current_key] = adata.obs.loc[adata.obs[current_key].isna(), parent_key]
    
    # Calculate normalized entropy for each new subcluster
    label_entropy = {}
    for leiden_value in adata.obs[current_key].unique():
        label_entropy[leiden_value] = cluster_entropy(adata, current_key, leiden_value, label_key)
    
    # Find clusters with high entropy for next level (using the normalized 0-1 threshold)
    high_entropy_clusters = [leiden for leiden, entropy in label_entropy.items() if entropy > label_entropy_thres]

    print(f"Found {len(high_entropy_clusters)} high-entropy clusters for next level")
    
    # Determine number of neighbors for next level
    next_n_neighbors = n_neighbors
    if level == 2:  # For level 3, use 15 neighbors
        next_n_neighbors = 15
    elif level >= 3:  # For levels 4 and 5, use 10 neighbors
        next_n_neighbors = 10
    
    next_resolution = resolution
    
    if level >= 3:
        next_resolution = 0.1
        
    # Recursively subcluster high-entropy clusters
    if level < 5 and high_entropy_clusters:
        recursive_subclustering(adata, high_entropy_clusters, level+1, next_n_neighbors, next_resolution, label_entropy_thres, current_key, label_key)
        
        
def top_mkgenes(adata, leiden, n):
    """
    Extracts the top N marker genes for each cluster based on previous rank_genes_groups calculation.
    Filters for positive log fold change and sorts by score.
    
    Parameters:
    -----------
    adata : AnnData
        The data object (RNA or ATAC gene scores).
    leiden : str
        The column in adata.obs containing cluster labels.
    -----------
    
    Returns:
    top_genes : Dictionary
        The dictionary containing top up-regulated genes for each cluster.
    """
    # Check if rank_genes_groups exists in uns
    if 'rank_genes_groups' not in adata.uns:
        raise ValueError("Please run sc.tl.rank_genes_groups() (or hf.calculate_t_stats) before running this function.")

    top_genes = {}
    
    # Ensure we iterate over string representations of clusters to match Scanpy's storage
    clusters = adata.obs[leiden].unique()
    
    for cluster in clusters:
        cluster_key = str(cluster)
        
        # Check if the cluster key exists in the analysis results
        if cluster_key in adata.uns['rank_genes_groups']['names'].dtype.names:
            
            # Extract data from structured array
            names = adata.uns['rank_genes_groups']['names'][cluster_key]
            logfoldchanges = adata.uns['rank_genes_groups']['logfoldchanges'][cluster_key]
            scores = adata.uns['rank_genes_groups']['scores'][cluster_key]
            pvals = adata.uns['rank_genes_groups']['pvals_adj'][cluster_key]
             
            # Create a DataFrame
            df = pd.DataFrame({'names': names, 'logfoldchanges': logfoldchanges, 'scores': scores, 'pvals': pvals})

            # Filter: logfoldchange > 0
            df = df[(df['logfoldchanges'] > 0) & (df['pvals'] < 0.05)]

            # Sort by scores descending and take top n
            top_df = df.sort_values(by='scores', ascending=False).head(min(n,df.shape[0]))

            # Store in dictionary
            top_genes[cluster] = top_df['names'].values

    # Print results
    for cluster, genes in top_genes.items():
        print(f'Cluster: {cluster}, Top {n} genes: {genes}')
        
    return top_genes