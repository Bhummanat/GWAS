# Run scripts in following sequence:
## Done with SLURM arrays
1. arrays.sh              ### Setup for SLURM arrays
2. vcf_subsetting.sh      ### Subsetting .vcf
3. subsetting_allele.R    ### Subsetting .vcf imports
4. joining_LGG_data.R     ### Joining SNPs dataset with covariates dataset
5. cox_modeling.R         ### Regularized coxph
## No SLURM arrays
6. FDR_and_KM.R           ### FDR-adjusted p-values and Kaplan Meier plot
7. roc_test.R             ### pROC (baseline vs SNP)
8. manhattan_qqplot.R     ### Combines statistics and produces Manhattan and qqplots
9. fine-mapping.Rmd
