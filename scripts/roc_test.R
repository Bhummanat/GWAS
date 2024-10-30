library(data.table)
library(tidyverse)
library(glmnet)
library(survival)
library(broom)
library(ggfortify)
library(qqman)
library(pROC)

##### Set the directory where the .rds files are stored ####
input_dir <- "/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/"

#### Import Files ####
LGG_allele_dat_sig <- readRDS(paste0(input_dir, 'LGG_allele_dat_sig.rds'))
LGG_pat <- readRDS(paste0(input_dir, 'LGG_pat.rds'))
gc()  # Run garbage collection after loading large files

dynamic_coxph_formula <- function(data, covariate_ext) { # data is each position; covariate_ext is extracted from each model
  covariates_all <- LGG_pat %>% colnames() %>% str_flatten(collapse = '|')
  covariates_matched <- str_extract(covariate_ext, covariates_all) %>% unique()
  
  covariates_to_include <- covariates_matched %>% 
    keep(~ n_distinct(data[[.x]], na.rm = TRUE) > 1)
  
  if ('GT' %in% covariate_ext && !'GT' %in% covariates_to_include) {
    return(NULL)
  }
  
  formula_text <- str_c('Surv(OS.time, OS) ~ ', str_c(covariates_to_include, collapse = '+')) %>% as.formula()
}

#### Generalized Model ####
LGG_pat_clean <- LGG_pat |> filter(!if_any(everything(), is.na))
x_imputed <- LGG_pat_clean %>% select(-c(OS.time, OS, GT)) %>% makeX()
y <- Surv(LGG_pat_clean$OS.time, LGG_pat_clean$OS)

cvfit <- cv.glmnet(x_imputed, y, family = "cox", alpha = 0.5, nfolds = 3)
lambda_1se <- cvfit$lambda.1se
coeff <- as.matrix(coef(cvfit, s = lambda_1se))
covariate_sig <- rownames(coeff)[coeff != 0]
covariate_generalized <- c('GT', covariate_sig)

rm(x_imputed, y, cvfit, coeff)  # Remove large objects not needed further
gc()  # Clear memory

#### ROC-test ####
run_coxph_per_position <- function(data, covariates_test, group_info) {
  position_value <- group_info$position
  data$position <- position_value
  formula <- dynamic_coxph_formula(data, covariates_test)
  
  if (is.null(formula)) return(NULL)  # Skip if formula is not valid
  
  model <- coxph(formula, data = data)
  used_data <- model.frame(model)
  excluded_rows <- data[!rownames(data) %in% rownames(used_data), ]
  
  risk_score <- predict(model, type = "risk")
  pat_included <- anti_join(data, excluded_rows, by = "bcr_patient_barcode")
  roc_pROC <- roc(pat_included$OS, risk_score)
  
  list(position = position_value, roc = roc_pROC)
}

#### AUC Calculations ####
coxph_baseline_sig <- LGG_allele_dat_sig %>%
  group_by(position) %>%
  group_map(~ run_coxph_per_position(.x, covariate_sig, .y))

rm(covariate_sig)  # Remove unneeded covariate variable
gc()  # Clear memory

coxph_generalized_sig <- LGG_allele_dat_sig %>%
  group_by(position) %>%
  group_map(~ run_coxph_per_position(.x, covariate_generalized, .y))

rm(covariate_generalized)  # Clear memory for final comparison
gc()  # Garbage collect again

roc_test_results <- map2(
  coxph_baseline_sig, 
  coxph_generalized_sig, 
  ~ {
    if (!is.null(.x) && !is.null(.y)) {
      baseline_roc <- .x$roc
      generalized_roc <- .y$roc
      roc.test(baseline_roc, generalized_roc, methods = 'bootstrap', boot.n = 1000)
    } else {
      NULL  # Skip if either result is NULL
    }
  }
)

p_values <- map_dbl(roc_test_results, ~ .x$p.value)

roc_test_summary <- tibble(
  position = map_chr(coxph_baseline_sig, ~ .x$position),
  p_value = p_values
)

write_csv(roc_test_summary, 'SNPs_AUC.csv')
