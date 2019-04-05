#!/bin/bash

# ABOUT: Script to perform a logistic regression model GWAS on 1000 Genomes Project (1KGP) phase 1 data using sex as trait.
# REQUIRED: 1000 Genomes Project phase 1 data (.bed, .bim, .fam) (see README)
# AUTHOR: Koen Rademaker
# DATE: 5 April 2019


###### Preparing the work environment and gathering data #####
# Set work directory
work_dir=#~/path/to/work_directory
# Set output directory
out_dir=#~/path/to/data_dir
# (Optional) Copy 1KGP data to work directory
cp 1kh_phase1_all* ${work_dir}
cd ${work_dir}


##### Preparing the data #####
# Recode data, randomly selecting 1% of variants
plink --bfile 1kg_phase1_all --thin 0.01 --recode --out PLINK_run
# Create .bed, .bim and .bim files for subsampled data
plink --bfile 1kg_phase1_all --thin 0.01 --make-bed --out PLINK_run
# Set sex as phenotype
awk '{print $1, $2, $3, $4, $5, $5}' PLINK_run.fam > tmp
mv tmp PLINK_run.fam


##### Running logistic regression model ######
# Run model with parameters MAF > 0.01 and exclusion of non-autosoal chromosomes
plink --bfile PLINK_run --logistic --maf 0.01 --not-chr 23-26 x y xy mt --out sex_autosomal_results


##### Formatting data for visualization #####
awk '{OFS="\t"; print $1, $2, $3, $9}' sex_autosomal_results.assoc.logistic > sex_autosomal_sum_stats.txt
