library(data.table)
library(tidyverse)
library(broom)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
chromosome <- args[1] # Format is chr#

# Joining
# SNPs Dataset
LGG_allele_filt <- readRDS(paste0('/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/', chromosome, '/LGG', chromosome, '_allele_filt.rds')) %>%
  as_tibble() %>%  # position, GT
  pivot_longer(contains('TCGA'), names_to = 'barcode', values_to = 'GT') %>% # Each position will have associated patients with GT values
  unite('position', `#CHROM`:POS) %>% 
  filter(!is.na(GT)) %>% # Per position, keep only patient that GT info is present
  mutate(barcode = str_sub(barcode, 1, 12)) %>% 
  group_by(position, barcode) %>%
  slice_head(n = 1) %>%  # Keep the first occurrence per patient-position pair (Some position took more than one sample)
  ungroup()

# Covariates datasets
patient_dat <- fread('/standard/vol169/cphg_ratan/ar7jq/TCGA/metadata/LGG/TCGA-CDR-SupplementalTableS1.tsv') %>% as_tibble() # age_at_initial_pathologic_diagnosis, gender, histological_grade, histological_type, OS, OS.time
patient_ancestry <- fread('/standard/vol169/cphg_ratan/ar7jq/TCGA/metadata/LGG/TCGAA_STRUCTURE.txt') %>% as_tibble() # race
patient_followup <- fread('/standard/cphg-RLscratch/syv3ap/rotation/filtered_TCGA_imports/clinical_PANCAN_patient_with_followup.csv.gz') %>% as_tibble() # TSS
aneuploidy <- fread('/standard/vol169/cphg_ratan/ar7jq/TCGA/metadata/LGG/ArmCallsAndAneuploidyScore_092817.txt') %>% as_tibble() %>% #aneuploidy, 1p/19q codeletion
  mutate(bcr_patient_barcode = str_sub(Sample, 1, 12)) %>% 
  group_by(bcr_patient_barcode) %>% 
  summarize(
    `Aneuploidy Score` = mean(`Aneuploidy Score`, na.rm = TRUE), # For each patient, `Aneuploidy Score` is average of their samples
    bcr_patient_barcode = first(bcr_patient_barcode),
    `1p` = first(`1p`),
    `19q` = first(`19q`)
  ) %>% 
  ungroup() %>% 
  mutate(
    percent_aneuploidy = `Aneuploidy Score`/39 * 100,
    `1p/19q_codeletion` = if_else((`1p` == -1 & `19q` == -1), 'present', 'absent')
  )
maf <- fread('/standard/vol169/cphg_ratan/ar7jq/TCGA/somatic/LGG/mc3.v0.2.8.PUBLIC.maf.gz') %>% as_tibble() %>% # somatic mutation count, IDH mutant
  select(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  mutate(
    bcr_patient_barcode = str_sub(Tumor_Sample_Barcode, 1, 12)) %>% 
  group_by(bcr_patient_barcode) %>%
  summarise(
    somatic_mut = n(),
    IDH_mut = if_else(sum(Hugo_Symbol %in% c('IDH1', 'IDH2')) > 0, 'mutant', 'wildtype'),
    TP53_mut = if_else(sum(Hugo_Symbol %in% c('TP53')) > 0, 'mutant', 'wildtype')
  )
cnv <- fread('/standard/cphg-RLscratch/syv3ap/rotation/filtered_TCGA_imports/TCGA_mastercalls.abs_segtabs.fixed.csv.gz') %>% as_tibble() %>% #Chr7 gain / Chr10 loss
  select(barcode, Chromosome, Modal_Total_CN) %>% 
  filter(Chromosome %in% c(7, 10)) %>% 
  group_by(barcode) %>% 
  summarize(
    `7G/10L` = if_else(any(Chromosome == 7 & Modal_Total_CN > 2) & any(Chromosome == 10 & Modal_Total_CN < 2), 
                       'present', 'absent')
  )
methylation <- fread('/standard/cphg-RLscratch/syv3ap/rotation/filtered_TCGA_imports/methylation_MGMT_promoter.csv') %>% as_tibble() %>% # MGMT promoter methylation
  summarize(across(contains('TCGA'), 
                   ~ mean(.x, na.rm = TRUE))) %>% 
  pivot_longer(everything(), names_to = 'barcode', values_to = 'avg_beta_val') %>%
  mutate(barcode = str_sub(barcode, 1, 12)) %>% 
  group_by(barcode) %>% 
  summarize(
    avg_beta_val = mean(avg_beta_val, na.rm = TRUE),
    barcode = first(barcode)
  ) %>% 
  ungroup() %>% 
  mutate(
    MGMT_prom_met = if_else(avg_beta_val > 0.3, 'methylated', 'unmethylated')
  )
tert_exp <- fread('/standard/vol169/cphg_ratan/ar7jq/TCGA/expression/LGG/TCGA.LGG.sampleMap%2FHiSeqV2_PANCAN.gz') %>% as_tibble() %>% 
  filter(sample == 'TERT') %>% 
  pivot_longer(contains('TCGA'), names_to = 'barcode', values_to = 'expression') %>% 
  mutate(
    barcode = str_sub(barcode, 1, 12),
  ) %>% 
  group_by(barcode) %>% 
  summarize(
    avg_expression = mean(expression, na.rm = TRUE), 
    barcode = first(barcode)
  ) %>% 
  ungroup() %>% 
  mutate(
    tert_shift_value = abs(min(avg_expression, na.rm = TRUE)),  # Calculate shift value
    tert_avg_expression = avg_expression + tert_shift_value      # Apply the shift
  )

patient <- left_join(patient_dat %>% select(bcr_patient_barcode, age_at_initial_pathologic_diagnosis, gender, histological_grade, histological_type, OS, OS.time), 
                     patient_ancestry %>% select('Patient ID', European, 'West African', 'Native American'), 
                     join_by(bcr_patient_barcode == 'Patient ID')) %>% 
  left_join(patient_followup %>% select(bcr_patient_barcode, tissue_source_site, tumor_location), 
            join_by(bcr_patient_barcode)) %>% 
  left_join(aneuploidy, join_by(bcr_patient_barcode)) %>% 
  left_join(maf, join_by(bcr_patient_barcode)) %>% 
  left_join(cnv, join_by(bcr_patient_barcode == barcode)) %>% 
  left_join(methylation, join_by(bcr_patient_barcode == barcode)) %>% 
  right_join(LGG_allele_filt, join_by(bcr_patient_barcode == barcode)) %>% # LGG_allele_filt will determine samples included
  left_join(tert_exp, join_by(bcr_patient_barcode == barcode)) %>% 
  filter(!is.na(OS.time)) %>% # Per position, keep only patient that GT info is present, keep only rows with valid OS.time
  group_by(position) %>% # In each position,
  filter(n_distinct(GT) > 1) %>% # there must be at least 2 genotype to compare
  group_by(position, GT) %>% # For each genotype in the position,
  filter(n_distinct(OS) > 1) %>% # there should be event variation (avoiding perfect separations and infinite coefficients)
  ungroup()%>%
  mutate(
    histological_grade = ordered(histological_grade, levels = c('G1', 'G2', 'G3', 'G4')), # Tumor-grade is ordered
    tissue_source_site = tissue_source_site %>%
      fct_other(keep = c("DU", "E1", "HT"), other_level = "Other") # Reducing data according to reference paper
  ) %>%
  mutate(across(                         # Change all other categoricals to factors
    c(GT, gender, histological_type, tumor_location, `7G/10L`, IDH_mut, TP53_mut, 
      `1p/19q_codeletion`, MGMT_prom_met),
    as_factor
  )) %>% 
  rename(                               
    'West_African' = `West African`,
    'Native_American' = `Native American`,
    'codeletion_1p19q' = `1p/19q_codeletion`,
    'chr_7G10L'= `7G/10L`
  )

saveRDS(patient, file = paste0('/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline/', chromosome, '/LGG', chromosome, '_allele_dat.rds'))

