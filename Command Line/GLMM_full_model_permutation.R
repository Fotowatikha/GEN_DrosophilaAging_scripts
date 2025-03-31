#!/usr/bin/env Rscript

# Author: Hamid Fotowatikha 
# Department: Laboratory of Genetics
# Date: 15-11-24
# Used for: Research Practice - Aging in Drosophila
# Title: Script to run the weighted GLMM on the permuted datasets

library(future)
library(furrr)
library(lmerTest)
library(lme4)
library(dplyr)
library(glmmTMB)
library(car)
library(future.apply)

#################### Chunk #################### 
# // Run the full  GLMM with beta distribution. Full model includes the age_group, and its interaction with fly_group (CE, CL)
# // This model runs on the permuted data
####################       ####################

#################### Chunk for COMMAND LINE #################### 
# // Run the FULL GLMM model with beta distribution 
# // Here we run it for the observed data and 20 permuted data frames
# // With multi threading
####################       ####################

# Set up parallel processing
options(future.globals.maxSize = 256000 * 1024^2)  # set MEM to 256GB
plan(multisession, workers = parallel::detectCores() - 32)  # Use available cores

# Define input and output directories
input_dir <- "/lustre/BIF/nobackup/fotow002/data_extreme_pheno/CRISP_OUT_variablePS_2025/GLM_GWAS/Permuted_data/Permuted_data_SoftWeight/"
output_dir <- "/lustre/BIF/nobackup/fotow002/data_extreme_pheno/CRISP_OUT_variablePS_2025/GLM_GWAS/GLM_output/GLM_results_permutation_SoftWeight/"
# Check if output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List all CSV files in the input directory
file_list <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# Function to process a single file
process_vcf_file <- function(file_path) {
  message("Processing file: ", file_path)

  # Read the data
  vcf_short <- read.csv(file_path, stringsAsFactors = FALSE)
  # Make the columns numerric
  vcf_short$weight <- as.numeric(vcf_short$weight)
  vcf_short$allele_frequency <- as.numeric(vcf_short$allele_frequency)

  # Ensure allele frequencies stay in the valid range (0 < y < 1)
  epsilon <- .Machine$double.eps
  vcf_short$allele_frequency <- ifelse(
    vcf_short$allele_frequency == 0, epsilon,
    ifelse(vcf_short$allele_frequency == 1, 1 - epsilon, vcf_short$allele_frequency)
  )
  rm(epsilon)

  # Remove NA values from allele frequency column
  vcf_short <- vcf_short[!is.na(vcf_short$allele_frequency), ]

  # Make nested lists for variants to be tested
  variant_list <- vcf_short %>%
    group_by(variant_id) %>%
    group_split() %>%
    future_map(identity)

  # Function to fit GLMM and extract p-values
  fit_glmm_extract_pvalues <- function(variant_data) {
    tryCatch({
      if (nlevels(as.factor(variant_data$fly_group)) < 2 || 
          nlevels(as.factor(variant_data$age_group)) < 2) {
        message("Skipping variant:", unique(variant_data$variant_id), " - Not enough factor levels")
        return(NULL)
      }
      
      # reduced GLMM model with beta distribution
      model <- glmmTMB(
        allele_frequency ~ fly_group + age_group*fly_group + (1 | biological_group) + (1 | biological_group:technical_group),
        data = variant_data, contrasts=list(age_group="contr.sum", fly_group="contr.sum"),
        family = beta_family(link = "logit"), 
        weights = variant_data$weight
      )
      
      # Extract p-values from Wald Chi-square Type3 test
      summary_model_Chi <- Anova(model, type = 3)
      p_value_fly_group_chi <- summary_model_Chi$`Pr(>Chisq)`[2] # Extract fly_group p-value fly group (CE vs CL)
      p_value_age_group_chi <- summary_model_Chi$`Pr(>Chisq)`[3] # Extract age_group p-value age group (young vs Old)
      p_value_age_group_fly_group_chi <- summary_model_Chi$`Pr(>Chisq)`[4]  # Extract age_group p-value interaction (Young vs Old in CE and CL)
      
      # Make required format of the dataframe and append values
      data.frame(
        CHROM = unique(variant_data$CHROM),
        POS = unique(variant_data$POS),
        variant_id = unique(variant_data$variant_id),
        REF = unique(variant_data$REF),
        ALT = unique(variant_data$ALT),
        VT = unique(variant_data$VT),
        FLANKSEQ = unique(variant_data$FLANKSEQ),
        NP = unique(variant_data$NP),
        VP = unique(variant_data$VP),
        DP = unique(variant_data$DP),
        CT = unique(variant_data$CT),
        VF = unique(variant_data$VF),
        EMstats = unique(variant_data$EMstats),
        p_value_fly_group = p_value_fly_group_chi,
        p_value_age_group = p_value_age_group_chi,
        p_value_age_fly_group = p_value_age_group_fly_group_chi, # interaction term
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("Error in variant:", unique(variant_data$variant_id))
      message(e)
      return(NULL)
    })
  }

  # Apply function in parallel
  test_results_list <- future_lapply(variant_list, fit_glmm_extract_pvalues, future.seed = TRUE)
  # Combine results into a single dataframe, filtering out failed models
  test_results <- do.call(rbind, test_results_list[!sapply(test_results_list, is.null)])
  # Replace NA p-values with 1
  test_results$p_value_fly_group[is.na(test_results$p_value_fly_group)] <- 1
  test_results$p_value_age_group[is.na(test_results$p_value_age_group)] <- 1
  test_results$p_value_age_fly_group[is.na(test_results$p_value_age_fly_group)] <- 1
  # Generate a unique filename based on the input file
  file_name <- basename(file_path)  # Extracts just the filename
  output_file <- file.path(output_dir, paste0("GLMM_results_", file_name))
  # Save results
  write.csv(test_results, output_file, row.names = FALSE)
  # Message fot Nohup
  message("Finished processing: ", file_name)
}

# Loop through all files in the directory, process them and run the GLLMM one by one
for (file in file_list) {
  process_vcf_file(file)
}

# Stop parallel processing
plan(sequential)
# Nohup message
cat("All files processed successfully without errors!!!\n")
