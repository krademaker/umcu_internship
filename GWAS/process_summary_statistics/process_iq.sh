#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process IQ (Intelligence Quotient) summary statistics data with (1) INFO > 0.9, (2) MAF > 0.01, (3) P < 0.05, (4) bi-allelic variants, (5) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Savage et al. (2018) IQ summary statistics, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 28 March 2019


# (1) File organization
summary_file="$TMPDIR"/iq/output_iq/summary.txt
output_file="$TMPDIR"/iq/output_iq/processed_iq_2018_sum_stats.txt
depict_file="$TMPDIR"/iq/output_iq/depict_iq_2018.txt
mkdir "$TMPDIR"/iq # Create work directory
mkdir "$TMPDIR"/iq/output_iq # Create output directory
cp $HOME/Koen/GWAS_data/control_gwas/iq_2018_sum_stats.txt "$TMPDIR"/iq # Copy summary statistics file to work directory
cp $HOME/1K/g1000_eur.synonyms "$TMPDIR"/iq # Copy synonyms file for 1KGP SNPs to work directory
cd "$TMPDIR"/iq # Move to work directory


# (2) Data processing
  # (2a) MAF filter
awk '$7>0.01 && $7<0.99' iq_2018_sum_stats.txt > iq_tmp_maf.txt # Filter out SNPs with MAF < 0.01
	# (2b) INFO filter
awk '$13>0.9' iq_tmp_maf.txt > iq_tmp_maf_info.txt # Filter out SNPs with INFO <= 0.9
	# (2c) P filter
awk '$11<0.05' iq_tmp_maf_info.txt > iq_tmp_maf_info_p.txt # Filter out SNPs with P >= 0.05
	# (2d) Bi-allelic filter
awk '$5=="g" || $5=="a" || $5=="c" ||$5=="t"' iq_tmp_maf_info_p.txt > iq_tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$6=="g" || $6=="a" || $6=="c" ||$6=="t"' iq_tmp_bialleles.txt > iq_tmp_maf_info_p_biallelic.txt # Filter out multi-allelic non-effect alleles
awk 'BEGIN { print "SNP\tUNIQUE_ID\tCHR\tPOS\tA1\tA2\tEAF_HRC\tZscore\tstdBeta\tSE\tP\tN_analyzed\tminINFO\tEffectDirection" } { print }' iq_tmp_maf_info_p_biallelic.txt > iq_tmp_maf_info_p_biallelic_header.txt # Re-add header to file
	# (2e) 1KGP filter
awk 'BEGIN { OFS=""; print "SNP" } FNR>2 { for (i=1; i<=NF; i++) print "rs",$i }' g1000_eur.synonyms > iq_tmp_rs_g1000_eur.synonyms # Reformat SNP synonyms to include 'rs' prefix
awk 'NR == FNR{c[$1]++;next};c[$1] > 0' iq_tmp_rs_g1000_eur.synonyms iq_tmp_maf_info_p_biallelic_header.txt > ${output_file} # Filter out SNPs not included in 1KGP list
	# (2f) DEPICT formatting
awk 'BEGIN { OFS = "\t" } FNR>1 {print $1,$3,$4,$11}' ${output_file} > iq_tmp_depict.txt # Extract columns 'SNP', 'CHR', 'BP' and 'P'
awk 'BEGIN { print "SNP\tChr\tPos\tP" }{ print }' iq_tmp_depict.txt > ${depict_file} # Rename columns to 'SNP', 'Chr', 'Pos', 'P' for final DEPICT file


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l iq_2018_sum_stats.txt >> ${summary_file} # Original file size
wc -l iq_tmp_maf.txt >> ${summary_file} # Effect of INFO filter
wc -l iq_tmp_maf_info.txt >> ${summary_file} # Effect of MAF filter
wc -l iq_tmp_maf_info_p.txt >> ${summary_file} # Effect of P filter
wc -l iq_tmp_maf_info_p_biallelic_header.txt >> ${summary_file} # Effect of bi-allelic filter
wc -l ${output_file} >> ${summary_file} # Processed file size
wc -l ${depict_file} >> ${summary_file} # DEPICT file size
	# (3b) Export output directory
cp -r "$TMPDIR"/iq/output_iq $HOME/Koen/GWAS_data
	# (3c) Clean temporary files
rm iq_tmp*
