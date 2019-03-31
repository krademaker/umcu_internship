#!/bin/bash

# ABOUT: Script for QC of GWAS summary statistics, generalized to handle various formats that data can be provided in
# REQUIRED: Summary statistics file, (optional) list of reference high-quality imputed SNPs, (optional) list of reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 31 March 2019


# (0) Decare input files and data columns
sum_stats=#/path/to/sum_stats
imputed_snps=#/path/to/imputed_snps
reference_snps=#/path/to/reference_snps

maf_col=# Set to MAF column
info_col=# Set to INFO column
info_snp_col=# Set to SNP column
p_col=# Set to P column
a1_col=# Set to A1 column
a2_col=# Set to A2 column

# (1) Perform QC(optionally write to temporary files)
awk -v col="$maf_col" '$col<0.01 && $col>0.99' ${sum_stats} # Filter MAF (default: >0.01)

awk -v col="$info_col" '$col>0.9' ${sum_stats} # Filter INFO score (default: >0.9)

awk -v col="$info_snp_col" '{print $col}' ${imputed_snps} > tmp_imputed_snps # Filter INFO score from reference list
mv tmp_imputed_snps ${imputed_snps}
awk 'NR == FNR {c[$1]++; next}; c[$1] > 0' ${imputed_snps} ${sum_stats}

awk -v col="$p_col" '$col<0.05' ${sum_stats} # Filter P-value (default: <0.05)

awk -v a1="$a1_col" -v a2="$a2_col" '($a1=="G" || $a1=="A" || $a1=="C" || $a1=="T") && ($a2=="G" || $a2=="A" || $a2=="C" || $a2=="T")' ${sum_stats} # Filter bi-allelic alleles
