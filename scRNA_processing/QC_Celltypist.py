#! python3
###STEPI+II: ######################################################
#1)The normal QC STEP
#2) Automatic cell type annotation using Celltypist
############################################################################################################
 
import pandas as pd
import scanpy as sc
import numpy as np
import os
import sys
from scipy.stats import median_abs_deviation
import scrublet as scr
from scipy.sparse import csr_matrix
import celltypist
from celltypist import models

tissue=sys.argv[1]
os.chdir(f"/storage/yangjianLab/chenchang/scRNA/{tissue}/raw_data")

#---------------------------------------------------------
print("Load data")
platform = sys.argv[2]

if platform=="10X":
    project = sys.argv[3] #PRJNA661032
    sample = sys.argv[4] #SAMN15967322
    path_filter = os.path.join(project, sample, 'outs/filtered_feature_bc_matrix')
    adata = sc.read_10x_mtx(path_filter)

elif platform=="10X_h5":
    project = sys.argv[3] #PRJNA661032
    sample = sys.argv[4] #SAMN15967322
    path_filter = os.path.join(project, f'raw_matrix/{project}_{sample}.h5')
    adata = sc.read_10x_h5(path_filter)
    adata.var_names_make_unique()

elif platform=="h5ad":
    project = sys.argv[3]
    sample = sys.argv[4]
    path_filter = os.path.join(project, f'raw_matrix/{project}_{sample}.h5ad')
    adata = sc.read_h5ad(path_filter)
    adata.raw=None
    adata.obs.drop(['orig.ident','nCount_RNA','nFeature_RNA'], axis=1, inplace=True)
    adata.var.drop('features', axis=1, inplace=True)

elif platform=="10X_multi":
    project = sys.argv[3] #PRJEB50594
    multi = sys.argv[4] #CRMulti
    sample = sys.argv[5] #DonorA
    path_filter = os.path.join(project, multi, "outs/per_sample_outs", sample, "count/sample_filtered_feature_bc_matrix")
    adata = sc.read_10x_mtx(path_filter, gex_only=True)
    sample = multi + "_" + sample

elif platform=="10X_csv":
    import anndata as ad
    project = sys.argv[3]
    sample = sys.argv[4]
    path_filter = f'{project}/raw_matrix/{project}_{sample}.counts.csv.gz'
    data_df = pd.read_csv(path_filter, index_col=0)
    sparse_matrix = csr_matrix(data_df)
    adata = ad.AnnData(X=sparse_matrix.transpose())
    adata.var_names = data_df.index
    adata.obs_names = data_df.columns

elif platform=="aggr":
    project = sys.argv[3] #EGAD00001007656
    sample = sys.argv[4] #LNA014_merged
    path_filter = os.path.join(project, sample, "outs/count/filtered_feature_bc_matrix")
    adata = sc.read_10x_mtx(path_filter)

elif platform=="multiaggr":
    project = sys.argv[3] #PRJEB45072
    sample = sys.argv[4] #P1
    path_filter = os.path.join(project, sample, "outs/count/filtered_feature_bc_matrix")
    adata = sc.read_10x_mtx(path_filter, gex_only=True)

elif platform=="freemuxlet":
    project = sys.argv[3] #PRJNA609655
    multi = sys.argv[4] #D1-D6_old
    sample = sys.argv[5] #0
    path_filter = os.path.join(project, multi, "outs/filtered_feature_bc_matrix")
    adata = sc.read_10x_mtx(path_filter)
    path_freemuxlet = os.path.join(project, multi, "freemuxlet", "barcode_sample_" + sample + ".csv")
    barcodes = pd.read_csv(path_freemuxlet, header=None)[0]
    adata = adata[barcodes, :]
    sample = multi + "_" + sample

elif platform=="singleron":
    project = sys.argv[3] #PRJNA661032
    sample = sys.argv[4] #SAMN15967322
    path_filter = os.path.join(project, sample, 'outs/filtered')
    adata = sc.read_10x_mtx(path_filter)

elif platform=="smartseq":
    import anndata as ad
    project = sys.argv[3]
    sample = sys.argv[4]
    path_filter = os.path.join(project, "raw_matrix", f"{project}_{sample}.genes.results")
    data_df = pd.read_csv(path_filter, sep="\t", index_col=0)
    sparse_matrix = csr_matrix(data_df)
    adata = ad.AnnData(X=sparse_matrix.transpose())
    adata.var_names = data_df.index
    adata.obs_names = data_df.columns    
    adata.var_names_make_unique()

elif platform=="BD":
    import anndata as ad
    project = sys.argv[3]
    sample = sys.argv[4]
    path_filter = f'{project}/raw_matrix/{project}_{sample}.counts.tsv.gz'
    data_df = pd.read_csv(path_filter, sep="\t", index_col=0)
    sparse_matrix = csr_matrix(data_df)
    adata = ad.AnnData(X=sparse_matrix.transpose())
    adata.var_names = data_df.index
    adata.obs_names = data_df.columns

elif platform=="dropseq":
    import anndata as ad
    project = sys.argv[3]
    sample = sys.argv[4]
    path_filter = f'{project}/raw_matrix/{project}_{sample}.dge.txt.gz'
    data_df = pd.read_csv(path_filter, sep="\t", index_col=0)
    sparse_matrix = csr_matrix(data_df)
    adata = ad.AnnData(X=sparse_matrix.transpose())
    adata.var_names = data_df.index
    adata.obs_names = data_df.columns

else:
    print("Please specify the correct platform.")

# Update column names
adata.obs_names = ["{}-{}-{}".format(name, project, sample) for name in adata.obs_names]
adata.obs["project"]=project
adata.obs["sample"]=sample

print("platform {}; project: {}; sample: {}".format(platform,project,sample))

#----------------------------------------------
# Basic filtering
print("Filtering")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
#adata.var["rb"] = adata.var_names.str.startswith(("RPS","RPL"))
#adata.var["hb"] = adata.var_names.str.startswith(("^HB[^(P)]"))
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=[20], log1p=True, inplace=True
)

def is_outlier(aadata, metric: str, nmads: int):
    M = aadata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
#    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 5) | (
    adata.obs["pct_counts_mt"] > 20
)
adata.obs.mt_outlier.value_counts()

print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

# Doublet detection using scrublet
print("Doublet detection")
ncell=adata.shape[0]
if ncell < 500:
    dbrate = 0.4 / 100
elif 500 <= ncell < 1000:
    dbrate = 0.4 / 100
elif 1000 <= ncell < 2000:
    dbrate = 0.8 / 100
elif 2000 <= ncell < 3000:
    dbrate = 1.6 / 100
elif 3000 <= ncell < 4000:
    dbrate = 2.3 / 100
elif 4000 <= ncell < 5000:
    dbrate = 3.1 / 100
elif 5000 <= ncell < 6000:
    dbrate = 3.9 / 100
elif 6000 <= ncell < 7000:
    dbrate = 4.6 / 100
elif 7000 <= ncell < 8000:
    dbrate = 5.4 / 100
elif 8000 <= ncell < 9000:
    dbrate = 6.1 / 100
elif 9000 <= ncell < 10000:
    dbrate = 6.9 / 100
elif ncell >= 10000:
    dbrate = 7.6 / 100

scrub = scr.Scrublet(adata.X,expected_doublet_rate=dbrate)
#doublet_scores, predicted_doublets = scrub.scrub_doublets()
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
adata.obs["doublet_score"]=doublet_scores
adata.obs["doublet_class"]=predicted_doublets
adata = adata[~adata.obs.doublet_class].copy()

ncell=adata.shape[0]
if ncell>=50:
    print(f'After filtering, {ncell} cells remaining.')
else:
    raise ValueError("Number of cells is less than 50, discard the sample.")

#----------------------------------------------
print("Automated cell type annotation using CellTypist")
adata_celltypist = adata.copy()  # make a copy of adata
sc.pp.normalize_total(adata_celltypist, target_sum=1e4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform
# make .X dense instead of sparse, for compatibility with celltypist:
adata_celltypist.X = adata_celltypist.X.toarray()
#models.download_models(model=["Human_Lung_Atlas"])
# please select the corresponding model for the tissue
model = models.Model.load(model="Human_Lung_Atlas.pkl")
predictions = celltypist.annotate(adata_celltypist, model=model, majority_voting=True)
predictions_adata = predictions.to_adata()
adata.obs["celltypist_pred"] = predictions_adata.obs.loc[adata.obs.index, "majority_voting"]
#adata.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[adata.obs.index, "conf_score"]

#----------------------------------------------
print("Writing out h5ad")
sample_h5ad_out = "/storage/yangjianLab/chenchang/scRNA/{}/h5ad/{}_{}_rawQC.h5ad".format(tissue,project, sample)
adata.write_h5ad(sample_h5ad_out)

print("Session finished")

