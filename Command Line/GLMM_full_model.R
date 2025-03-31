#!/usr/bin/env Rscript

# Author: Hamid Fotowatikha 
# Department: Laboratory of Genetics
# Date: 15-11-24
# Used for: Research Practice - Aging in Drosophila
# Title: Script to run the weighted GLMM on the full dataset post pre-processing of CRISP output (SoftWeight_vcf_short_full.csv)

#################### Chunk #################### 
# // Apply the weighted GLMM model in on SoftWeight_vcf_short_full.csv
# // With multi threading
# // This script will outputGLMM_results_full_model_SoftWeight.csv
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
# // Run the full  GLMM with beta distribution. Full model includes the age_group and its interaction with fly_group (CE, CL)
# // Correction for multiple testing is not included here, as we this will be done in the next part (see the followup script)
####################       ####################

# Set up parallel processing
options(future.globals.maxSize = 256000 * 1024^2) # set mem on 256Gb
plan(multisession, workers = parallel::detectCores() - 64)  # Use all but one core

# Path to vcf_short_full.csv
vcf_short <- read.csv("/lustre/BIF/nobackup/fotow002/data_extreme_pheno/CRISP_OUT_variablePS_2025/GLM_GWAS/vcf_short_FROM_R/SoftWeight_vcf_short_full.csv")

# Shift 0s to a small value and 1s to just below 1 to make it compatible with the GLMM model and Beta distribution
# .Machine$double.eps is the smallest number so that 1 + eps > 1 in double-precision floars in R
# - If allele frequency is exactly 0, replace it with epsilon (smallest positive distinguishable number)
# - If allele frequency is exactly 1, replace it with (1 - epsilon) to avoid numerical issues
# - All other values remain unchanged in allele_frequencies
epsilon <- .Machine$double.eps
vcf_short$allele_frequency <- ifelse(
  vcf_short$allele_frequency == 0, epsilon,
  ifelse(vcf_short$allele_frequency == 1, 1-epsilon, vcf_short$allele_frequency)
)

# Remove all rows with allele_frequency columns equal to NA, representing ultra low_depth regions for some sample
# Note that this line is no longer in use, since we applied CRISP methodology filtering (CT, MAC48 adn DP360 filtering)
vcf_short <- vcf_short[!is.na(vcf_short$allele_frequency), ]

# make nested lists forvairants to be tested
variant_list <- vcf_short %>% # use vcf_short1 to test small dataset
  group_by(variant_id) %>%
  group_split() %>%
  future_map(identity)

# Function to fit GLMM and extract p-values
fit_glmm_extract_pvalues <- function(variant_data) {
  tryCatch({
    # Make sure that fly_group and age_group have at least 2 levels
    # Since we remove NA rows in the allele frequencies column in vcf_short, they will also remove the fly_group rows, so this line is usless and only used for testing
    if (nlevels(as.factor(variant_data$fly_group)) < 2 || 
        nlevels(as.factor(variant_data$age_group)) < 2) {
      message("Skipping variant:", unique(variant_data$variant_id), " - Not enough factor levels")
      return(NULL)  # Skip this variant
    }

    # Fit the model GLMM model with a beta distribution
    model <- glmmTMB(
      allele_frequency ~ fly_group + age_group*fly_group + (1 | biological_group) + (1 | biological_group:technical_group),
      data = variant_data, contrasts=list(age_group="contr.sum", fly_group="contr.sum"),
      family = beta_family(link = "logit"), 
      weights = variant_data$weight
    )

    # Extract p-values from Wald z-test
    #summary_model <- as.data.frame(summary(model)$coefficients$cond)
    #p_value_fly_group_z <- summary_model$`Pr(>|z|)`[2]  
    #p_value_age_group_z <- summary_model$`Pr(>|z|)`[3]  
    #p_value_age_group_fly_group_z <- summary_model$`Pr(>|z|)`[4]  
    
    # Extract p-values from Wald Chi-square Type3 test
    summary_model_Chi <- Anova(model, type = 3)
    p_value_fly_group_chi <- summary_model_Chi$`Pr(>Chisq)`[2] # Extract fly_group p-value fly group (CE vs CL)
    p_value_age_group_chi <- summary_model_Chi$`Pr(>Chisq)`[3] # Extract age_group p-value age group (young vs Old)
    p_value_age_group_fly_group_chi <- summary_model_Chi$`Pr(>Chisq)`[4]  # Extract age_group p-value interaction (Young vs Old in CE and CL)

    # Make a dataframe for storing results from columns we want
    data.frame(
      CHROM = unique(variant_data$CHROM),
      POS = unique(variant_data$POS),
      variant_id = unique(variant_data$variant_id),
      REF = unique(variant_data$REF),
      ALT = unique(variant_data$ALT),
      AF = unique(variant_data$AF),
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
  }, error = function(e) { # handle those fucking errors! (this is no longer a fucking problem in the GLMM with beta distribution)
    message("Error in variant:", unique(variant_data$variant_id)) # show the error in the given variant
    message(e)
    return(NULL)
  })
}

# Apply function in parallel
test_results_list <- future_lapply(variant_list, fit_glmm_extract_pvalues, future.seed = TRUE)
# Combine results into a single dataframe, filtering out failed models
test_results <- do.call(rbind, test_results_list[!sapply(test_results_list, is.null)])
# Change the NAN to 1 in case of failed test
test_results$p_value_fly_group[is.na(test_results$p_value_fly_group)] <- 1
test_results$p_value_age_group[is.na(test_results$p_value_age_group)] <- 1
test_results$p_value_age_fly_group[is.na(test_results$p_value_age_fly_group)] <- 1
# Apply multiple testing correction
#test_results$adjusted_p_value_fly_group <- p.adjust(test_results$p_value_fly_group, method = "bonferroni")
#test_results$adjusted_p_value_age_group <- p.adjust(test_results$p_value_age_group, method = "bonferroni")
#test_results$adjusted_p_value_age_fly_group <- p.adjust(test_results$p_value_age_fly_group, method = "bonferroni")
# Stop parallel processing
plan(sequential)
# Save or view the results
write.csv(test_results, "/lustre/BIF/nobackup/fotow002/data_extreme_pheno/CRISP_OUT_variablePS_2025/GLM_GWAS/GLM_output/GLMM_results_full_model_SoftWeight.csv", row.names = FALSE)

# Nohup messege
cat("Everything ran succesfully without errors!!!.\n")
