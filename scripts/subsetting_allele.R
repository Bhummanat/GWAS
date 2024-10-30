library(data.table)
library(tidyverse)
library(broom)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
chromosome <- args[1]

# Load the chromosome-specific VCF file
vcf_file <- paste0('/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/', chromosome, '_LGG_filtered.vcf.gz')
vcf1 <- fread(vcf_file, skip = 'CHROM') %>% as_tibble()

# Subsetting
vcf1 <- vcf1 %>% 
  mutate(
    info_split = str_split(INFO, ';'), # Since number of pieces varies
    AF = map_chr(info_split, ~ .x[1]), # Split and select index as follows
    AC = map_chr(info_split, ~ .x[length(.x) - 1]),
    AN = map_chr(info_split, ~ .x[length(.x)])
  )

vcf1 <- vcf1 %>% 
  mutate(
    AF = gsub(".*=", "", AF), # Clean-up to isolate numerical strings
    AC = gsub(".*=", "", AC), 
    AN = gsub(".*=", "", AN)
  )

vcf1$AF <- as.double(vcf1$AF) # Strings to integer
vcf1$AC <- as.integer(vcf1$AC) 
vcf1$AN <- as.integer(vcf1$AN)

vcf_allele_filt <- vcf1 %>% 
  filter(AN >= 200) %>% # How many alleles sequenced
  filter(AF > 0.05 & AF < 0.95) # Control for allele variance

vcf_allele_filt <- vcf_allele_filt %>%  
  rename_with(
    .fn = ~ str_extract(., "TCGA-[^_]*-[^_]*"),  # Extract only the first "TCGA-.*-.*" part
    .cols = contains("TCGA")  # Apply only to columns that have "TCGA" in their name
  ) %>% 
  mutate(
    `#CHROM` = `#CHROM`,
    POS = POS,
    across( # Selecting columns and changing #|# to # format
      .cols = contains("TCGA"), # Homo-reference_0, hetero_1, homo-alt_2
      .fns = ~ case_when(
        str_extract(., '[01][|/][01]') %in% c('0/0', '0|0') ~ 0,
        str_extract(., '[01][|/][01]') %in% c('0/1', '0|1', '1/0', '1|0') ~ 1,
        str_extract(., '[01][|/][01]') %in% c('1/1', '1|1') ~ 2,
        TRUE ~ NA_real_  # Catch any unexpected or missing values
      )
    ),
    .keep = "used"
  )

output_file <- paste0('/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/LGG', chromosome, '_allele_filt.rds')
saveRDS(vcf_allele_filt, file = output_file)


