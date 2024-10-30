library(data.table)
library(tidyverse)
library(broom)
library(survival)
library(ggfortify)

# Set the directory where the .rds files are stored
input_dir <- "/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/"

# Get a list of all the .rds files for each chromosome
LGG_allele_dat_parts <- list.files(input_dir, pattern = 'LGGchr.*_allele_dat.rds', full.names = TRUE, recursive = TRUE)
snp_basic_parts <- list.files(input_dir, pattern = "chr.*_snp_basic.csv.gz", full.names = TRUE, recursive = TRUE)
snp_generalized_parts <- list.files(input_dir, pattern = "chr.*_snp_generalized.rds", full.names = TRUE, recursive = TRUE)

#### Global FDR Correction for Significant SNPs ####

# Step 1: Apply FDR correction for basic SNPs and save significant ones
snp_basic <- map_dfr(snp_basic_parts, fread) %>%
  mutate(p.adjust = p.adjust(p.value, method = 'BH')) %>%
  filter(p.adjust < 0.1)
write_delim(snp_basic, file = 'snp_sig_basic.csv.gz', delim = ',')
sig_snps_basic <- snp_basic$snp_position

# Step 2: Apply FDR correction for generalized SNPs and save significant ones
snp_generalized <- map_dfr(snp_generalized_parts, readRDS) %>%
  mutate(p.adjust = p.adjust(p.value, method = 'BH')) %>%
  filter(p.adjust < 0.1)
write_delim(snp_generalized, file = 'snp_sig_generalized.csv.gz', delim = ',')
sig_snps_generalized <- snp_generalized$snp_position

#### Filter LGG Allele Data for Significant SNPs ####

# Load allele data in manageable batches and filter for significant SNPs
LGG_allele_dat_sig <- map_dfr(LGG_allele_dat_parts, function(file) {
  readRDS(file) %>%
    filter(position %in% sig_snps_generalized)
})
saveRDS(LGG_allele_dat_sig, file = 'LGG_allele_dat_sig.rds')
gc()

##### Process Patient Data for ROC Test ####
LGG_pat <- LGG_allele_dat_sig %>%
  filter(OS.time > 0) %>%
  group_by(bcr_patient_barcode) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-c(bcr_patient_barcode, `Aneuploidy Score`, `1p`, `19q`, avg_beta_val, position, avg_expression, tert_shift_value))
saveRDS(LGG_pat, file = 'LGG_pat.rds')
rm(LGG_pat)
gc()

#### Generate KM Plots in Batches ####
# Ensure KM_Plots directory exists
dir.create("/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/KM_Plots", recursive = TRUE, showWarnings = FALSE)

# Lower batch size to manage memory
batch_size <- 2  

position_batches <- LGG_allele_dat_sig %>%
  group_by(position) %>%
  group_split() %>%
  split(ceiling(seq_along(.) / batch_size))

map(position_batches, ~ {
  walk(.x, ~ {
    if (nrow(.x) == 0) return(NULL)
    
    km_snps <- survfit(Surv(OS.time, OS) ~ GT, data = .x)
    plot <- autoplot(km_snps, conf.int = TRUE) + 
      labs(title = str_c('Position: ', unique(.x$position)))
    
    ggsave(filename = str_c('KM_Plot_', unique(.x$position), '.png'), 
           plot = plot, 
           path = '/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/KM_Plots',
           dpi = 100)
    
    rm(plot)
    gc()
  })
  gc()
})
