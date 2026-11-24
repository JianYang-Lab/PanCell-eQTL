"""Create raw cell-by-gene expression matrices and donor QC tables for one tissue.

Use raw counts in X and donor IDs in obs["project_sample"] matching the ancestry
FAM files. Group cells by celltypist_final and Ancestry_pred. Write each eligible
group to the tissue's h5ad_celltype directory.
"""

import argparse
from pathlib import Path
from urllib.parse import quote

import anndata as ad
import numpy as np
import pandas as pd


def prepare_annotations(adata, tissue_dir):
    tissue = tissue_dir.name
    obs = adata.obs
    
    if "project_sample" not in obs or ("Ancestry_pred" not in obs and "sample" in obs):
        sample = obs["project"].astype(str) + "_" + obs["sample"].astype(str)
        obs["project_sample"] = pd.Categorical(sample)

    if "Ancestry_pred" not in obs:
        path = tissue_dir / "sc-eQTL_R" / "covariates" / f"{tissue}_all_sample_cov.txt"
        ancestry = pd.read_csv(path, sep="\t", dtype=str).set_index("sample")["Ancestry_pred"]
        obs["Ancestry_pred"] = pd.Categorical(obs["project_sample"].astype(str).map(ancestry))

    if "celltypist_final" not in obs:
        paths = sorted((tissue_dir / "h5ad_merge" / "cell_type_final").glob("*.txt"))
        if not paths:
            raise ValueError(
                "Supply celltypist_final in the H5AD or barcode labels in h5ad_merge/cell_type_final/"
            )
        labels = pd.concat(pd.read_csv(p, sep="\t", index_col=0).iloc[:, 0] for p in paths)
        obs["celltypist_final"] = pd.Categorical(labels.reindex(obs.index))


def write_matrix(adata, obs, prefix):
    counts = obs["project_sample"].value_counts()
    donors = counts.index[counts >= 5]
    if len(donors) < 30:
        return

    cells = obs.index[obs["project_sample"].isin(donors)].to_numpy()
    groups = adata.obs.iloc[cells].groupby("project_sample", observed=True)
    n_expressing = np.zeros(adata.n_vars, dtype=int)
    for rows in groups.indices.values():
        n_positive = np.asarray((adata.X[cells[rows]] > 0).sum(axis=0)).ravel()
        n_expressing += n_positive >= 5

    # Keep genes expressed in at least five cells in at least 10% of eligible donors.
    keep = 10 * n_expressing >= len(donors)
    adata[cells, keep].write_h5ad(f"{prefix}.h5ad")
    qc = groups[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].mean().round(4)
    qc.columns = ["avg_total_counts", "avg_n_genes", "avg_mt_pct"]
    qc["cell_count"] = groups.size()
    qc.to_csv(f"{prefix}_cov.csv")
    print(f"{prefix.name}: {len(donors)} donors, {keep.sum()} genes")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tissue", type=Path, help="Tissue directory, e.g. Blood")
    parser.add_argument(
        "--input", type=Path, help="Raw-count H5AD (default: TISSUE/h5ad_merge/TISSUE_merged_new.h5ad)"
    )
    args = parser.parse_args()
    tissue_dir = args.tissue.resolve()
    tissue = tissue_dir.name
    adata = ad.read_h5ad(args.input or tissue_dir / "h5ad_merge" / f"{tissue}_merged_new.h5ad")
    prepare_annotations(adata, tissue_dir)
    out = tissue_dir / "h5ad_celltype"
    out.mkdir(exist_ok=True)

    labels = adata.obs[["Ancestry_pred", "celltypist_final", "project_sample"]].reset_index(drop=True)
    for ancestry, obs in labels.groupby("Ancestry_pred", observed=True):
        fam = tissue_dir / "sc-eQTL_R" / "genotype" / f"{tissue}_{ancestry}_indivpruned.fam"
        if not fam.is_file():
            print(f"Skipping {ancestry}: missing {fam.name}")
            continue
        genotyped = pd.read_csv(fam, sep=r"\s+", header=None, usecols=[1], dtype=str)[1]
        eligible = obs[
            obs["project_sample"].astype(str).isin(genotyped) & (obs["celltypist_final"] != "Discard")
        ]
        for cell_type, cells in eligible.groupby("celltypist_final", observed=True):
            label = quote(cell_type, safe=" ")  # Cell type labels may contain a slash.
            write_matrix(adata, cells, out / f"{tissue}_rawQC_{label}_{ancestry}")


if __name__ == "__main__":
    main()
