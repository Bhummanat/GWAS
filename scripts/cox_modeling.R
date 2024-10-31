.libPaths("~/R/goolf/4.3") # Or your preferred custom library directory
library(data.table)
library(tidyverse)
library(glmnet)
library(survival)
library(broom)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
chromosome <- args[1] # Format is chr#

#### Import Dataset ####
LGG_allele_dat <- readRDS(paste0('/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/', chromosome, '/LGG', chromosome, '_allele_dat.rds')) %>% 
  mutate(GT = as.integer(as.character(GT)))

##### For feature selection, we are not considering SNPs, therefore we reduce dataset to consider only patient-relevants ####
LGG_pat <- LGG_allele_dat %>%
  filter(OS.time > 0) %>% # 'OS.time == 0' does not work with model
  group_by(bcr_patient_barcode) %>% # Extract one datapoint of each patient (patient datapoint was duplicated for SNPs)
  slice_head(n = 1) %>% 
  ungroup() |> 
  select(-c(bcr_patient_barcode, `Aneuploidy Score`, `1p`, `19q`,  # Reducing by removing intermediate variables
            avg_beta_val, position, avg_expression, tert_shift_value))

#### Dynamic cox formula ####
dynamic_coxph_formula <- function(data, covariate_ext){ # data is each position; covariate_ext is extracted from each model
  
  # Significant covariates (String Wrangling)
  covariates_all <- LGG_pat %>% colnames() %>% str_flatten(collapse = '|')
  covariates_matched <- str_extract(covariate_ext, covariates_all) %>%  # covariate_ext is defined by the basic or generalized model
    unique()
  
  # Filters
  covariates_to_include <- covariates_matched %>% 
    keep(~ n_distinct(data[[.x]], na.rm = TRUE) > 1) # This removed variables that are invalid for the current data
  
  # Check if 'GT' is in covariate_ext; Meaning entered model depends on GT
  if ('GT' %in% covariate_ext) {
    # Don't run model in current position if 'GT' is not included; This is required as there might not be variation of genotype in the current data
    if (!'GT' %in% covariates_to_include) {
      return(NULL) # NULL output to flag this; Should be handled by whatever used this 'dynamic_coxph_formula' function
    }
  }
  
  # Formula with valid covariates
  formula_text <- str_c('Surv(OS.time, OS) ~ ', str_c(covariates_to_include, collapse = '+')) %>% 
    as.formula()
}

#### For basic model, all covariates are considered ####
covariate_full <- LGG_pat %>% select(-c(OS, OS.time)) %>% colnames() # All covariates assigned for basic coxph model

# Basic coxph over for all snps
coxph_basic <- LGG_allele_dat %>% # Basic coxph model is created from entire patient dataset
  group_by(position) %>% # For each snps
  group_map(
    ~ {
      formula <- dynamic_coxph_formula(.x, covariate_full) # Creating formula based on datapoints in current position with all covariates
      
      # If GT is not included and NULL flag returned
      if (is.null(formula)){
        return(tibble(term = 'GT not included', snp_position = unique(.x$position))) # Skip the rest (No coxph for current position)
      }
      
      model <- coxph(formula, data = .x, control = coxph.control(iter.max = 10000)) # Model using dynamic formula per position
      
      tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% 
        mutate(snp_position = unique(.x$position)) %>% # Attach position to output
        filter(str_detect(term, 'GT')) # Choose data associated with GT
    },
    .keep = TRUE
  )

# Extracting significant SNPs from basic modeling
snp_basic <- bind_rows(coxph_basic) %>% 
  filter(!is.na(p.value))

#### Generalized Model ####
# Setup for cross-validation
LGG_pat_clean <- LGG_pat |> filter(!if_any(everything(), is.na)) # Removed data with NA as NA doesn't seem to have correlation
x_imputed <- LGG_pat_clean %>% select(-c(OS.time, OS, GT)) %>% makeX() # makeX(na.impute = TRUE) can fill NA with average as an alternative
y <- Surv(LGG_pat_clean$OS.time, LGG_pat_clean$OS)

# Extracting significant covariates (non-zero coefficients at lambda.1se)
cvfit <- cv.glmnet(x_imputed, y, family = "cox", alpha = 0.5, nfolds = 3) # Elastic-net
lambda_1se <- cvfit$lambda.1se
coeff <- as.matrix(coef(cvfit, s = lambda_1se))
covariate_sig <- rownames(coeff)[coeff != 0]
covariate_generalized <- c('GT', covariate_sig)
write(covariate_generalized, '/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/covariate_generalized.txt')

# Generalized coxph over for all snps
coxph_generalized <- LGG_allele_dat %>% # Generalized coxph model is created from entire patient dataset
  group_by(position) %>% # For each snps
  group_map(
    ~ {
      formula <- dynamic_coxph_formula(.x, covariate_generalized) # Creating formula based on datapoints in current position with significant covariates
      
      # If GT is not included and NULL flag returned
      if (is.null(formula)){
        return(tibble(term = 'GT not included', snp_position = unique(.x$position))) # Skip the rest (No coxph for current position)
      }
      
      model <- coxph(formula, data = .x, control = coxph.control(iter.max = 10000)) # Model using dynamic formula per position
      
      tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% 
        mutate(snp_position = unique(.x$position)) %>% # Attach position to output
        filter(str_detect(term, 'GT')) # Choose data associated with GT
    },
    .keep = TRUE
  )

# Extracting significant SNPs from generalized modeling
snp_generalized <- bind_rows(coxph_generalized) %>% 
  filter(!is.na(p.value))

# Define the base output directory for the chromosome
output_dir <- paste0('/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/', chromosome, '/')

# Write to file: SNPs from basic model
write_delim(snp_basic, file = paste0(output_dir, chromosome, '_snp_basic.csv.gz'), delim = ',')

# Write to file: SNPs from generalized model
write_delim(snp_generalized, file = paste0(output_dir, chromosome, '_snp_generalized.csv.gz'), delim = ',')
saveRDS(snp_generalized, file = paste0(output_dir, chromosome, '_snp_generalized.rds'))


