# GWAS Pipeline in R

A pipeline for conducting genome-wide association studies (GWAS) using unix & R.

## Table of Contents
- [Introduction](#introduction)
- [Parts](#Parts)
- [Input Data Required](#InputDataRequired)

## Introduction
Pipeline for performing genome-wide association studies (GWAS) using R. It supports data preprocessing, statistical analysis, and visualization of results.

## Parts
- Subsetting and filtering imputed germline vcf
- SNP filtering based on Allele sequenced and frequency
- Merging with covariates datasets
- Elastic net feature selection & coxph
- Visualization of FDR-significant SNPs with kaplan-meier curve
- roc-test of FDR-significant SNPs
- Visualization of results with Manhattan and QQ plots

## InputDataRequired
- .vcf file containing population of interest
- .info file containing imputation information
- list of patients in 0_TCGA-.*-.*_TCGA-.*-.* line-delimited format

## Self-Reminder
Setup config.R where environment variables is set (So others don't need to use my directories)
