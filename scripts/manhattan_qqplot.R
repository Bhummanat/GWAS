library(data.table)
library(tidyverse)
library(broom)
library(qqman)

# Set the directory where the .rds files are stored
input_dir <- "/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/"

# Get a list of all the .rds files for each chromosome
snp_generalized_parts <- list.files(input_dir, pattern = "chr.*_snp_generalized.rds", full.names = TRUE, recursive = TRUE) # CHECK RECURSIVE

# Read all SNP data files and combine them
snp_generalized <- map_dfr(snp_generalized_parts, readRDS)

# Load ROC-test Data
roc_test <- fread('/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/SNPs_AUC.csv') %>%
  rename_with(~ str_replace(., 'p_value', 'p-roc'))

# Joining Datasets (Includes: CHR, BP, snp = 'snp' or AUC-sig, p-roc, conf. ind, exponentiate)
gwas_res <- left_join(snp_generalized, roc_test, join_by(snp_position == position)) %>%
  separate(snp_position, into = c('CHR', 'BP'), sep = '_')

# Write to file gwas_res (csv.gz and rds)
saveRDS(gwas_res, 'gwas_res.rds')
write_delim(gwas_res, 'gwas_res.csv.gz', delim = ',')

# Manhattan Plot
gwas_res <- readRDS('/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/gwas_res.rds') # Troubleshoot

png("/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/GWAS_manhattan_plot_rs34988193.png", width = 3000, height = 1200, res = 150)
gwas_res %>% 
  filter(CHR == '15') %>%  # Optional chromosome filtering
  select(CHR, BP, p.value) %>%
  mutate(
    CHR = case_match(
      CHR,
      'X' ~ '23',
      'Y' ~ '24',
      .default = CHR),
    p.value = p.value + 1e-10, # Set upper y_lim to 10
    CHR = as.integer(CHR),
    BP = as.integer(BP),
    SNP = as.character(BP) # Don't need this in data if SNP assigned to each position prior
  ) %>%
  manhattan(p = 'p.value', main = "Manhattan plot of GWAS p-values", cex = 0.4, annotatePval = 1e-5, annotateTop = FALSE,
            highlight = snpsOfInterest) # Optional highlights (Must enter snps of interest in console)
dev.off()

# QQplot
png("GWAS_qqplot.png", width = 3000, height = 1200, res = 150)
qq(gwas_res$p.value, main = "Q-Q plot of GWAS p-values")
dev.off()
