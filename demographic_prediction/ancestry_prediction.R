#!/usr/bin/env Rscript

# ==============================================================================
# Name: 5-ancestry_pred.R
# Description: This script uses Principal Components (PCs) from a combined
#              genotype dataset (sample data + 1KGP reference) to predict
#              the ancestry of the samples. It includes optional batch
#              correction using Harmony and trains a multinomial logistic
#              regression model for prediction.
# ==============================================================================

# --- Load Libraries ---
suppressMessages({
    library(tidyverse)
    library(caret)
    library(nnet)
    library(harmony)
    library(argparse)
})

# --- Function Definitions ---

#' Load and preprocess PCA and demographic data.
#' @param bed_dir Directory containing PCA eigenval/eigenvec files.
#' @param demo_file Path to the sample demographics file.
#' @param ref_pop_file Path to the 1KGP reference population file.
#' @return A processed data frame with PCs and population labels.
load_and_process_data <- function(bed_dir, demo_file, ref_pop_file) {
    cat("Step 1: Loading and processing data...\n")

    # Load PCA results
    eigenval <- sqrt(read.table(file.path(bed_dir, "merged_ref_pruned_pca.eigenval"))[, 1])
    eigenvec <- read.table(file.path(bed_dir, "merged_ref_pruned_pca.eigenvec"))
    
    # Scale eigenvectors
    num_pcs <- ncol(eigenvec) - 2
    for (k in 1:num_pcs) {
        eigenvec[, 2 + k] <- eigenvec[, 2 + k] * eigenval[k]
    }
    vec <- as.data.frame(eigenvec[, -1])
    colnames(vec) <- c("sample", paste0("PC", 1:num_pcs))

    # Load and merge demographic data
    pop_table <- read_delim(demo_file, delim = '\t', show_col_types = FALSE)
    ref_pop_table <- read.table(ref_pop_file, col.names = c("sample", "ancestry"))
    
    all_pop <- bind_rows(
        pop_table %>% select(sample, ancestry),
        ref_pop_table
    ) %>% distinct(sample, .keep_all = TRUE)

    # Combine PCA with population data
    vec_with_pop <- vec %>%
        left_join(all_pop, by = "sample") %>%
        mutate(
            REF = ifelse(sample %in% ref_pop_table$sample, "Reference", "Sample"),
            POP = ifelse(ancestry %in% c("Hispanic or Latino", "Hispanic/Latino", "Jewish", "Middle Eastern", "Other"), NA, ancestry)
        )

    cat("...Data loading complete.\n\n")
    return(vec_with_pop)
}

#' Run Harmony for batch correction between sample and reference PCs.
#' @param pca_data The data frame containing the original PCs.
#' @return A data frame with batch-corrected PCs.
run_harmony_correction <- function(pca_data) {
    cat("Step 2: Running Harmony for batch correction...\n")
    
    pc_matrix <- pca_data %>% select(starts_with("PC"))
    
    harmony_corrected <- HarmonyMatrix(
        data_mat = pc_matrix,
        meta_data = pca_data,
        vars_use = "REF",
        do_pca = FALSE # PCs are already calculated
    )
    
    corrected_df <- as.data.frame(harmony_corrected)
    # Copy over metadata
    corrected_df$sample <- pca_data$sample
    corrected_df$REF <- pca_data$REF
    corrected_df$POP <- pca_data$POP
    
    cat("...Harmony correction complete.\n\n")
    return(corrected_df)
}

#' Perform 5-fold cross-validation of the ancestry prediction model.
#' @param pc_data The data frame with PCs to use for validation.
run_cross_validation <- function(pc_data) {
    cat("Step 3: Performing 5-fold cross-validation...\n")
    set.seed(100) # for reproducibility

    # Prepare data for cross-validation
    cv_data <- pc_data %>% filter(!is.na(POP))
    
    # Define the population levels for consistent factor encoding
    pop_levels <- c("EUR", "EAS", "AFR", "AMR", "SAS", "Other")
    cv_data$POP <- factor(cv_data$POP, levels = pop_levels)
    
    folds <- createFolds(cv_data$POP, k = 5, list = TRUE, returnTrain = FALSE)
    accuracies <- numeric(5)

    for (i in 1:5) {
        test_indices <- folds[[i]]
        train_data <- cv_data[-test_indices, ]
        test_data <- cv_data[test_indices, ]

        model <- multinom(
            POP ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + REF,
            data = train_data,
            decay = 5e-4,
            maxit = 200,
            trace = FALSE
        )

        predictions <- predict(model, test_data)
        
        # Ensure predictions and actual values have the same factor levels
        predictions_factor <- factor(predictions, levels = pop_levels)
        actual_factor <- factor(test_data$POP, levels = pop_levels)
        
        cm <- confusionMatrix(predictions_factor, actual_factor)
        accuracies[i] <- cm$overall["Accuracy"]
        cat("  Fold", i, "accuracy:", accuracies[i], "\n")
    }

    mean_accuracy <- mean(accuracies)
    sd_accuracy <- sd(accuracies)
    cat("  Mean five-fold accuracy:", mean_accuracy, "\n")
    cat("  Standard deviation for accuracy:", sd_accuracy, "\n")
    cat("...Cross-validation complete.\n\n")
}


#' Train a model and predict ancestry for all samples.
#' @param pc_data The data frame with PCs to use for training/prediction.
#' @return A data frame with sample IDs and their predicted ancestries.
predict_ancestry <- function(pc_data) {
    cat("Step 4: Training final model and predicting ancestry for all samples...\n")

    # Prepare training data (known ancestries, mostly from reference)
    train_data <- pc_data %>% filter(!is.na(POP))

    # Train a multinomial logistic regression model on the first 10 PCs + REF status
    model <- multinom(
        POP ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + REF,
        data = train_data,
        decay = 5e-4,
        maxit = 200,
        trace = FALSE
    )

    # Predict on all samples
    predictions <- predict(model, newdata = pc_data)
    probabilities <- predict(model, newdata = pc_data, type = "probs")
    
    # Get the max probability for confidence scoring
    max_probs <- apply(as.data.frame(probabilities), 1, max)

    results <- data.frame(
        sample = pc_data$sample,
        ancestry_pred = as.vector(predictions),
        ancestry_prob = max_probs
    )
    
    cat("...Prediction complete.\n\n")
    return(results)
}

#' Save the final ancestry predictions.
#' @param predictions The data frame with prediction results.
#' @param original_demographics The original sample demographics table.
#' @param output_dir The directory to save the output files.
save_results <- function(predictions, original_demographics, output_dir) {
    cat("Step 5: Saving results...\n")
    
    # Filter for only the project's samples
    project_samples <- predictions %>%
        filter(sample %in% original_demographics$sample)

    # Set low-confidence predictions to NA
    project_samples$ancestry_pred[project_samples$ancestry_prob < 0.5] <- NA
    
    # Combine with original demographics
    final_output <- original_demographics %>%
        left_join(project_samples %>% select(sample, ancestry_pred, ancestry_prob), by = "sample")

    # Define output paths
    pred_only_path <- file.path(output_dir, paste0(tissue, "_sample_ancestry_predicted.txt"))
    full_demo_path <- file.path(output_dir, paste0(tissue, "_sample_demographics_with_ancestry.txt"))

    # Write files
    write.table(
        final_output %>% select(sample, ancestry_pred, ancestry_prob),
        file = pred_only_path,
        row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE
    )
    write.table(
        final_output,
        file = full_demo_path,
        row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE
    )
    
    cat(paste("  Predicted ancestries saved to:", pred_only_path, "\n"))
    cat(paste("  Full demographics with predictions saved to:", full_demo_path, "\n"))
    cat("...Results saved.\n")
}


# --- Main Execution ---

# Setup command line argument parser
parser <- ArgumentParser(description = "Predict sample ancestry using PCA data.")
parser$add_argument("--tissue", type = "character", required = TRUE,
                    help = "Name of the tissue (e.g., Colon).")
# TODO: Update these default paths to match your project structure
parser$add_argument("--bed_dir", type = "character", 
                    default = "/path/to/your/project/{tissue}/VCF_TopMedimputed/filtered",
                    help = "Directory containing the merged and pruned PCA results.")
parser$add_argument("--demo_file", type = "character",
                    default = "/path/to/your/project/{tissue}/demographics/sample_demographics.txt",
                    help = "Path to the sample demographics file.")
parser$add_argument("--ref_pop_file", type = "character",
                    default = "/path/to/your/reference_panel/1KGP_3202_sample_ancestry.txt",
                    help = "Path to the 1000 Genomes reference population file.")
parser$add_argument("--output_dir", type = "character",
                    default = "./ancestry_results",
                    help = "Directory to save the output files.")

args <- parser$parse_args()

# Substitute {tissue} placeholder in paths
tissue <- args$tissue
bed_dir <- gsub("\\{tissue\\}", tissue, args$bed_dir)
demo_file <- gsub("\\{tissue\\}", tissue, args$demo_file)

# Create output directory
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Workflow ---
original_demographics <- read_delim(demo_file, delim = '\t', show_col_types = FALSE)

pca_data <- load_and_process_data(bed_dir, demo_file, args$ref_pop_file)
#corrected_pcs <- run_harmony_correction(pca_data)
run_cross_validation(pca_data) # New cross-validation step
predictions <- predict_ancestry(pca_data)
save_results(predictions, original_demographics, args$output_dir)

cat("\n--- Ancestry prediction workflow completed successfully! ---\n")
