#!/bin/bash
#SBATCH -n 1
#SBATCH -t 15:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to annotate summary statistics with rs IDs (reference SNPs) when only SNP positions are given
# REQUIRED: UCSC SNP table (see README for details), GWAS summary statistics file
# AUTHOR: Koen Rademaker
# DATE: 28 March 2019


# (1) File organization
output_file="$TMPDIR"/rs_id_annotation/output_rs_id_annotation/annotated_scz_2018_sum_stats.txt
mkdir "$TMPDIR"/rs_id_annotation # Create work directory
mkdir "$TMPDIR"/rs_id_annotation/output_rs_id_annotation # Create output directory
cp "$HOME"/Koen/GWAS_data/scz_gwas/scz_2018_sum_stats.txt "$TMPDIR"/rs_id_annotation # Copy SCZ summary statistics file to work directory
cp $HOME/Koen/UCSC_SNP_annotation/ucsc_snp_table.txt "$TMPDIR"/rs_id_annotation # Copy UCSC SNP annotation file to work directory
cd "$TMPDIR"/rs_id_annotation # Move to work directory

sum_stats="$TMPDIR"/rs_id_annotation/scz_2018_sum_stats.txt
ucsc_snps="$TMPDIR"/rs_id_annotation/ucsc_snp_table.txt
restructured_ucsc_snps="$TMPDIR"/rs_id_annotation/tmp_ucsc.txt
annotated_sum_stats="$TMPDIR"/rs_id_annotation/output_rs_id_annotation/rs_id_annotation_scz_2018_sum_stats.txt


# (2) Restructure UCSC SNP table
awk 'BEGIN { OFS=""; print "SNP\tRS" } $2==$3 {gsub("chr", ""); print $1, ":", $2, "\t", $4}' ${ucsc_snps} > ${restructured_ucsc_snps} # Extract items with 1bp length, reformat SNP name

# (3) Annotate rs IDs for SNPs
  # (3a) Create comparison file for cut and uncut SNP names
awk -v n=2 'BEGIN { OFS="\t"; print "SNP_cut\tSNP_uncut" } FNR>1 { s = gensub("^(([^:]*:){"n"}).*$", "\\1", 1); sub(".$","",s); print s,$1}' ${sum_stats} > rs_tmp_snp_comparison.txt
  # (3b) Cut allele information from SNP column
cut -f1,2 -d':' ${sum_stats}.txt | awk '{print $1}' > rs_tmp_cut_snps.txt
  # (3c) Annotate rs IDs to SNPs
awk 'BEGIN { OFS="\t" } NR==FNR{A[$1];next}$1 in A' rs_tmp_cut_snps.txt ${restructured_ucsc_snps} > rs_tmp_cut_snps_rs.txt
  # (3d) Restore uncup SNPs
awk 'BEGIN { OFS="\t" } FNR==NR{a[$1]=$2;next}{print $0,a[$1]?a[$1]:"NA"}' rs_tmp_cut_snps_rs.txt rs_tmp_snp_comparison.txt > rs_tmp_combined_snps_rs.txt
  # (3e) Export uncut SNPs with rs ID annotation
awk 'BEGIN { OFS="\t", print "SNP\tRS" } $3 != "NA" {print $2, $3}' rs_tmp_combined_snps_rs.txt > ${annotated_sum_stats}


# (4) File export and cleaning
	# (4a) Export output directory
cp -r "$TMPDIR"/rs_id_annotation/output_rs_id_annotation $HOME/Koen/UCSC_SNP_annotation
	# (4b) Clean temporary files
rm rs_tmp*
