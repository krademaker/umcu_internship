#!/bin/bash

# ABOUT: Script to annotate summary statistics with rsIDs (reference SNPs) when only SNP positions are given
# REQUIRED: UCSC SNP table (see README for details), GWAS summary statistics file
# AUTHOR: Koen Rademaker
# DATE: 27 March 2019

# (0) Set up variables
sum_stats=#~/path/to/sum_stats.txt
ucsc_snps=#~/path/to/ucsc_snp_table.txt
restructured_ucsc_snps=#~/path/to/tmp_ucsc.txt
annotated_sum_stats=#~/path/to/annotated_sum_stats.txt

# (1) Restructure UCSC SNP table
awk 'BEGIN { OFS=""; print "SNP\tRS" } $2==$3 {gsub("chr", ""); print $1, ":", $2, "\t", $4}' ${ucsc_snps} > ${restructured_ucsc_snps} # Extract items with 1bp length, reformat SNP name

# (2) Annotate rs IDs for SNPs
awk -v n=2 'BEGIN { OFS="\t"; print "SNP_cut\tSNP_uncut" } FNR>1 { s = gensub("^(([^:]*:){"n"}).*$", "\\1", 1); sub(".$","",s); print s,$1}' ${sum_stats} > tmp_snp_comparison.txt  # Create comparison file for cut and uncut SNP names
cut -f1,2 -d':' ${sum_stats}.txt | awk '{print $1}' > tmp_cut_snps.txt # Cut allele information from SNP column
awk 'BEGIN { OFS="\t" } NR==FNR{A[$1];next}$1 in A' tmp_cut_snps.txt ${restructured_ucsc_snps} > tmp_cut_snps_rs.txt # Annotate rs IDs to SNPs
awk 'BEGIN { OFS="\t" } FNR==NR{a[$1]=$2;next}{print $0,a[$1]?a[$1]:"NA"}' tmp_cut_snps_rs.txt tmp_snp_comparison.txt > tmp_combined_snps_rs.txt # Restore uncup SNPs
awk 'BEGIN { OFS="\t", print "SNP\tRS" } $3 != "NA" {print $2, $3}' tmp_combined_snps_rs.txt > ${annotated_sum_stats} # Export uncut SNPs with rs ID annotation

# (3) Clean temporary files
rm *tmp*
