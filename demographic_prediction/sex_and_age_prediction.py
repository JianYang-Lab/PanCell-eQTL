#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Predicts sample demographics (sex and age) from single-cell gene expression data.

This script performs the following steps:
1.  Loads and preprocesses an AnnData object and demographic metadata.
2.  Performs cross-validation to find optimal hyperparameters for a Random
    Forest model, then predicts sex using sex-linked DEGs.
3.  Performs cross-validation for the age prediction model and reports performance,
    then predicts age using an AutoGluon TabularPredictor trained on age-related DEGs.
4.  Saves the predicted sex and age for all samples to specified output files.
"""

import os
import argparse
import warnings
import shutil
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.model_selection import train_test_split, GridSearchCV, KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
from autogluon.tabular import TabularDataset, TabularPredictor

# Suppress common warnings for a cleaner output
warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None


def load_and_preprocess_data(h5ad_path, demographics_path, exclude_path, rename_path):
    """Loads and preprocesses the core datasets."""
    print("Step 1: Loading and preprocessing data...")

    # Load AnnData and metadata
    adata = sc.read_h5ad(h5ad_path)
    meta_df = pd.read_csv(demographics_path, sep="\t")
    meta_df.set_index("sample", inplace=True)

    # Exclude specified barcodes (e.g., aneuploid cells)
    if exclude_path and os.path.exists(exclude_path):
        bc_exclude = pd.read_csv(exclude_path, header=None)[0].values
        adata = adata[~adata.obs_names.isin(bc_exclude)].copy()
        print(f"  Excluded {len(bc_exclude)} barcodes.")

    # Rename samples to match VCF IDs if a mapping file is provided
    if rename_path and os.path.exists(rename_path):
        rename_map = pd.read_csv(rename_path, sep=r'\s+')
        adata.obs['project_sample'] = adata.obs['project_sample'].astype(str)
        rename_dict = dict(zip(rename_map['h5ad'], rename_map['vcf']))
        adata.obs['project_sample'] = adata.obs['project_sample'].map(rename_dict).fillna(adata.obs['project_sample'])
        print(f"  Renamed samples based on '{os.path.basename(rename_path)}'.")

    # Standardize gender labels in metadata
    meta_df['gender'] = meta_df['gender'].astype(str).str.lower()
    meta_df.loc[meta_df['gender'].isin(["f", "female"]), 'gender'] = "female"
    meta_df.loc[meta_df['gender'].isin(["m", "male"]), 'gender'] = "male"
    
    print("...Data loading complete.\n")
    return adata, meta_df


def create_pseudobulk(adata, group_by_key='project_sample'):
    """Creates a pseudobulk expression matrix by averaging expression per sample."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    expression_df = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        index=adata.obs.index,
        columns=adata.var_names
    )
    
    pseudobulk_X = expression_df.groupby(adata.obs[group_by_key]).mean()
    return pseudobulk_X

def run_sex_prediction_cv(X_train, y_train):
    """Performs 5-fold cross-validation to find the best hyperparameters."""
    print("  Running 5-fold cross-validation for hyperparameter tuning...")
    clf = RandomForestClassifier(random_state=42)
    
    param_grid = {
        'n_estimators': [150, 300],
        'max_depth': [None, 10, 20],
    }
    
    cv = GridSearchCV(clf, param_grid, cv=5, scoring='accuracy', n_jobs=-1)
    cv.fit(X_train, y_train)
    
    print(f"  Best parameters found: {cv.best_params_}")
    print(f"  Best cross-validated accuracy: {cv.best_score_:.4f}")
    
    return cv.best_params_

def predict_sex(adata, meta_df, sex_deg_path, output_dir):
    """Predicts sex using a RandomForest model and saves the results."""
    print("Step 2: Predicting sex...")
    
    # Load sex DEGs and filter data
    sex_deg_df = pd.read_table(sex_deg_path)
    sex_genes = sex_deg_df.loc[sex_deg_df.chr == "chrX", "HUGO_gene_id"].values.tolist() + ["RPS4Y1"]
    
    valid_genes = [gene for gene in sex_genes if gene in adata.var_names]
    adata_sex = adata[:, valid_genes].copy()

    # Create pseudobulk expression matrix for sex genes
    X_sex_pseudobulk = create_pseudobulk(adata_sex)
    
    # Merge with metadata to get labels
    X_meta_sex = pd.merge(X_sex_pseudobulk, meta_df['gender'], left_index=True, right_index=True, how='left')
    X_meta_sex['gender'].fillna("NA", inplace=True)

    # Prepare data for ML model
    train_df = X_meta_sex[X_meta_sex['gender'].isin(['female', 'male'])]
    X_train_full = train_df.drop('gender', axis=1)
    y_train_full = train_df['gender']
    
    X_inference = X_meta_sex.drop('gender', axis=1)

    # Perform cross-validation to find best parameters
    best_params = run_sex_prediction_cv(X_train_full, y_train_full)

    # Train final model on all available labeled data with the best parameters
    print("  Training final Random Forest model with best parameters...")
    final_clf = RandomForestClassifier(random_state=42, **best_params)
    final_clf.fit(X_train_full, y_train_full)
    
    # Predict sex for all samples
    predicted_genders = final_clf.predict(X_inference)
    gender_predictions_s = pd.Series(predicted_genders, index=X_inference.index, name='gender_pred')

    # Save predictions
    output_path = os.path.join(output_dir, "predicted_sex.txt")
    gender_predictions_s.to_csv(output_path, sep='\t', index_label='sample')
    print(f"  Sex predictions saved to '{output_path}'.")
    
    # Return the original metadata DataFrame with the new predictions merged in
    meta_df_with_sex = pd.merge(meta_df, gender_predictions_s, left_index=True, right_index=True, how='left')
    
    print("...Sex prediction complete.\n")
    return meta_df_with_sex

def run_age_prediction_cv(train_data, output_dir):
    """Performs 5-fold cross-validation for the AutoGluon age model."""
    print("  Running 5-fold cross-validation for age prediction model...")
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    mae_scores = []

    for i, (train_index, val_index) in enumerate(kf.split(train_data)):
        print(f"    Processing Fold {i+1}/5...")
        train_fold_data = train_data.iloc[train_index]
        val_fold_data = train_data.iloc[val_index]
        
        fold_model_path = os.path.join(output_dir, f"autogluon_age_model_cv_fold_{i+1}")

        predictor = TabularPredictor(
            label='age',
            problem_type='regression',
            eval_metric='mean_absolute_error',
            path=fold_model_path
        ).fit(
            train_data=train_fold_data,
            presets='medium_quality', # Use a faster preset for CV
            num_cpus=10,
            time_limit=600 # 10 minute time limit per fold
        )
        
        performance = predictor.evaluate(val_fold_data)
        mae = performance['mean_absolute_error'] * -1 # AutoGluon makes it negative
        mae_scores.append(mae)
        print(f"    Fold {i+1} Mean Absolute Error: {mae:.4f}")
        
        # Clean up the fold-specific model directory to save space
        shutil.rmtree(fold_model_path, ignore_errors=True)

    mean_mae = np.mean(mae_scores)
    std_mae = np.std(mae_scores)
    
    print(f"\n  Average cross-validated Mean Absolute Error: {mean_mae:.4f} (+/- {std_mae:.4f})\n")


def predict_age(adata, meta_df_with_sex, age_deg_paths, output_dir):
    """Predicts age using AutoGluon and saves the results."""
    print("Step 3: Predicting age...")

    # Load age DEGs from multiple files
    age_genes = set()
    for f in age_deg_paths:
        if 'agingatlas' in f:
            age_genes.update(pd.read_csv(f)['Symbol'].tolist())
        elif 'genage' in f:
            age_genes.update(pd.read_csv(f)['symbol'].tolist())
    
    valid_age_genes = list(age_genes.intersection(adata.var_names))
    adata_age = adata[:, valid_age_genes].copy()

    # Create pseudobulk matrix for age genes
    X_age_pseudobulk = create_pseudobulk(adata_age)

    # Merge with metadata (age and predicted gender)
    X_meta_age = pd.merge(X_age_pseudobulk, meta_df_with_sex[['age', 'gender_pred']], left_index=True, right_index=True, how='left')

    # Prepare training and inference data
    age_known_mask = pd.to_numeric(X_meta_age['age'], errors="coerce").notna()
    train_df = X_meta_age[age_known_mask].copy()
    train_df['age'] = train_df['age'].astype(float)
    
    inference_df = X_meta_age.drop('age', axis=1)

    # Perform cross-validation on the training data
    run_age_prediction_cv(train_df, output_dir)

    # Train final AutoGluon model on all available training data
    train_data = TabularDataset(train_df)
    inference_data = TabularDataset(inference_df)

    model_path = os.path.join(output_dir, "autogluon_age_model_final")
    print(f"  Training final AutoGluon model (models saved to '{model_path}')...")
    predictor = TabularPredictor(
        label='age',
        problem_type='regression',
        eval_metric='mean_absolute_error',
        path=model_path
    ).fit(
        train_data=train_data,
        presets='good_quality',
        num_cpus=10,
        time_limit=3600 # Set a time limit of 1 hour for the final model
    )

    # Predict age for all samples
    age_predictions = predictor.predict(inference_data)
    age_predictions.name = "predicted_age"
    
    # Combine known ages with predicted ages for a final output file
    final_age_df = train_df[['age']].copy()
    final_age_df = pd.merge(final_age_df, age_predictions, left_index=True, right_index=True, how='right')
    final_age_df['final_age'] = final_age_df['age'].fillna(final_age_df['predicted_age'])

    output_path = os.path.join(output_dir, "predicted_age.txt")
    final_age_df['final_age'].to_csv(output_path, sep='\t', index_label='sample')
    print(f"  Final age predictions saved to '{output_path}'.")
    print("...Age prediction complete.\n")


def main():
    """Main function to orchestrate the prediction workflow."""
    parser = argparse.ArgumentParser(description="Predict demographics from gene expression data.")
    parser.add_argument('--h5ad_file', required=True, help='Path to the input AnnData (.h5ad) file.')
    parser.add_argument('--demographics_file', required=True, help='Path to the sample demographics TSV file.')
    parser.add_argument('--sex_deg_file', required=True, help='Path to the sex DEG list file.')
    parser.add_argument('--age_deg_files', required=True, nargs='+', help='Paths to one or more age DEG list files.')
    parser.add_argument('--output_dir', required=True, help='Directory to save output files and models.')
    parser.add_argument('--exclude_barcodes', default=None, help='(Optional) Path to a text file of barcodes to exclude.')
    parser.add_argument('--rename_file', default=None, help='(Optional) Path to a file for renaming samples.')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # --- Workflow ---
    adata, meta_df = load_and_preprocess_data(
        args.h5ad_file, args.demographics_file, args.exclude_barcodes, args.rename_file
    )
    
    meta_with_sex = predict_sex(
        adata, meta_df, args.sex_deg_file, args.output_dir
    )
    
    predict_age(
        adata, meta_with_sex, args.age_deg_files, args.output_dir
    )
    
    print("Workflow finished successfully!")


if __name__ == "__main__":
    main()
