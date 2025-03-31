#!/usr/bin/env Rscript

# Author: Hamid Fotowatikha 
# Department: Laboratory of Genetics
# Date: 15-11-24
# Used for: Research Practice - Aging in Drosophila
# Title: Script to make permuted datasets from the "SoftWeight_vcf_short_full_for_permutation.csv" after BH FDR correction from the first GLMM run

#################### Chunk #################### 
# // Make permutations on vcf_short
# // With multi threading
####################       ####################

library(future)
library(furrr)
library(lmerTest)
library(lme4)
library(dplyr)
library(glmmTMB)
library(car)
library(future.apply)

#################### Chunk #################### 
# // Make 20 Permutation according to Hoesjed et al., 2019: Jha et al. (2015)
# // This is required to correct for multiple testing and determine a significance threshold for differentiated SNPs
# // Pseudo-randomization is carried out under two criteria: (1) a sample could not retain its original label, and (2) replicates of a selection regime could not be reassigned as replicates of another regime.  
# // Here, the original vcf_short.csv will be used for permutations, followed by recalculation of weights for each variant within each permuted data
# // This is a more efficient code without per line loops
####################       ####################

# Set multithreading
options(future.globals.maxSize = 256000 * 1024^2)  # 256GB memory limit
plan(multisession, workers = parallel::detectCores() - 64)  # Use half of available cores

# Load original csv data or only significant ones to permute
vcf_short_permute <- read.csv("/lustre/BIF/nobackup/fotow002/data_extreme_pheno/CRISP_OUT_variablePS_2025/GLM_GWAS/vcf_short_FROM_R/SoftWeight_vcf_short_full_for_permutation.csv")

# Define output directory
output_dir <- "/lustre/BIF/nobackup/fotow002/data_extreme_pheno/CRISP_OUT_variablePS_2025/GLM_GWAS/Permuted_data/Permuted_data_SoftWeight/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to permute groups without retaining the same (fly_group, age_group)
permute_one_variant <- function(variant_data) {
  unique_combinations <- unique(variant_data[, c("fly_group", "age_group", "biological_group", "technical_group")])
  
  # Shuffle while ensuring no row retains the same (fly_group, age_group)
  # This is more ffcient than per line looping to check for unique options
  repeat {
    shuffled_combinations <- unique_combinations[sample(nrow(unique_combinations)), ]
    # Make sure no row keeps the same (fly_group, age_group)
    mismatch <- !(
      (unique_combinations$fly_group == shuffled_combinations$fly_group) & 
      (unique_combinations$age_group == shuffled_combinations$age_group)
    )
    if (all(mismatch)) break  # If all are valid, then stop repeating
  }
  # Assign the new values
  variant_data[, c("fly_group", "age_group", "biological_group", "technical_group")] <- shuffled_combinations
  
  return(variant_data)
}

# Function to calculate GLMMM weights
calculate_weight <- function(data) {
  data %>%
    group_by(variant_id) %>%
    mutate(
      weight = {
        short_living <- mean(allele_frequency[fly_group == "short_living"], na.rm = TRUE)
        long_living <- mean(allele_frequency[fly_group == "long_living"], na.rm = TRUE)
        allele_diff <- abs(short_living - long_living)
        max_frequency <- pmin(short_living, long_living)
        (allele_diff) * (1 - max_frequency)
      }
    ) %>%
    ungroup() # we ungroup in case of downstream dyplr:: funtion
}

# Number of permutations
num_permutations <- 20
# Main Permutation Loop
for (i in 1:num_permutations) {
  cat(sprintf("Generating permutation %d...\n", i))
  
  # Unique seed per iteration
  set.seed(Sys.time() %>% as.numeric() + i)
  
  # Copy original dataset
  permuted_df <- vcf_short_permute
  # Parallel processing of permutation
  permuted_df <- future_map_dfr(
    split(permuted_df, permuted_df$variant_id), 
    permute_one_variant, 
    .progress = TRUE
  )
  # Calculate weights (single-threaded)
  permuted_df <- calculate_weight(permuted_df)
  
  # Save permutation as CSV
  file_name <- paste0("permutation_", i, ".csv")
  file_path <- file.path(output_dir, file_name)
  write.csv(permuted_df, file = file_path, row.names = FALSE)
  
  cat(sprintf("Saved permutation %d to %s\n", i, file_path))
}

# Reset CPUs
plan(sequential)
cat("All permuted datasets have been saved successfully.\n")
