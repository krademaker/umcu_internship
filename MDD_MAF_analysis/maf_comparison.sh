#!/bin/bash

# ABOUT: Extract allele frequencies for both case and control MDD (Major Depressive Disorder) populations for further analysis
# AUTHOR: Koen Rademaker
# DATE: 25 March 2019


# Download summary statistics data (Wray et al. 2018)
wget https://www.med.unc.edu/pgc/results-and-downloads/data-use-agreement-forms/copy2_of_mdd2018_ex23andme%20_data_download_agreement/resolveuid/a9451724b03741a4bf708c0014d7ff12 --no-check-certificate -O mdd_2018_sum_stats.txt.gz
# Unzip data
gunzip mdd_2018_sum_stats.txt.gz
# Extract columns with case and control allele frequencies
awk 'BEGIN { OFS="\t" } FNR>1 { print $6, $7 }' mdd_2018_sum_stats.txt > tmp_maf.txt
# Add custom header to file
awk 'BEGIN { print "MAF_cases\tMAF_controls" }{ print }' tmp_maf.txt > maf_comparison.txt
# Remove temporary file
rm tmp_maf.txt
