#!/bin/bash

# ABOUT: Script for QC of GWAS summary statistics, generalized to handle various formats that data can be provided in
# REQUIRED: Summary statistics file, (optional) list of reference high-quality imputed SNPs, (optional) list of reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 31 March 2019


# DECLARE INPUT FILES AND DATA COLUMNS
sum_stats=#/path/to/sum_stats
imputed_snps=#/path/to/imputed_snps
reference_snps=#/path/to/reference_snps

maf_col=# Set to MAF column
info_col=# Set to INFO column
info_snp_col=# Set to SNP column
p_col=# Set to P column
a1_col=# Set to A1 column
a2_col=# Set to A2 column
snp_col=# Set to SNP (rsID) column
chr_col=# Set to chromosome column
pos_col=# Set to (base-pair) position column


# PERFORM QC STEPS (OPTIONAL: WRITE TO TEMPORARY FILES)
awk -v col="$maf_col" '$col<0.01 && $col>0.99' ${sum_stats} # Filter MAF (default: >0.01)

awk -v col="$info_col" '$col>0.9' ${sum_stats} # Filter INFO score (default: >0.9)

awk -v col="$info_snp_col" '{print $col}' ${imputed_snps} > tmp_imputed_snps # Filter INFO score from reference list
mv tmp_imputed_snps ${imputed_snps}
awk 'NR == FNR {c[$1]++; next}; c[$1] > 0' ${imputed_snps} ${sum_stats}

awk -v col="$p_col" '$col<0.05' ${sum_stats} # Filter P-value (default: <0.05)

awk -v a1="$a1_col" -v a2="$a2_col" '($a1=="G" || $a1=="A" || $a1=="C" || $a1=="T") && ($a2=="G" || $a2=="A" || $a2=="C" || $a2=="T")' ${sum_stats} # Filter bi-allelic alleles

grep 'rs' ${sum_stats} # Select SNPs with rsID

awk -v n=1 '{s = gensub("^(([^:]*:){"n"}).*$", "\\1", 1); sub(".$","",s); print s}' ${sum_stats} # Extract rsID from full SNP name

awk 'OFS=""; FNR>2 { for (i=1; i<=NF; i++) print "rs",$i }' ${reference_snps} > tmp_reference_snps # Filter reference SNPs (default: filtered MAGMA SNP synonyms, see README)
mv tmp_reference_snps ${reference_snps}
awk -v col="$snp_col" 'NR == FNR{c[$1]++;next};c[$col] > 0' ${reference_snps} ${sum_stats}

awk -v snp="$snp_col" -v chr="$chr_col" -v pos="$pos_col" -v p="$p_col" 'BEGIN {OFS="\t"; print "SNP\tChr\tPos\tP"} FNR>1 {print $snp,$chr,$pos,$p}' ${sum_stats} # Convert to DEPICT input format
