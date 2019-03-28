#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process BMI (Body Mass Index) summary statistics data with (1) INFO > 0.9, (2) MAF > 0.01, (3) P < 0.05, (4) bi-allelic variants, (5) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Yengo et al. (2018) BMI summary statistics, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATA: 26 March 2019


# (1) File organization
summary_file="$TMPDIR"/bmi/output_bmi/summary.txt
output_file="$TMPDIR"/bmi/output_bmi/processed_bmi_2018_sum_stats.txt
depict_file="$TMPDIR"/bmi/output_bmi/depict_bmi_2018.txt
mkdir "$TMPDIR"/bmi # Create work directory
mkdir "$TMPDIR"/bmi/output_bmi # Create output directory
cp $HOME/Koen/GWAS_data/control_gwas/bmi_2018_sum_stats.txt "$TMPDIR"/bmi # Copy summary statistics file to work directory
cp $HOME/1K/g1000_eur.synonyms "$TMPDIR"/bmi # Copy synonyms file for 1KGP SNPs to work directory
cd "$TMPDIR"/bmi # Move to work directory


# (2) Data processing
	# (2a) MAF filter
awk '$6>0.01 && $6<0.99' bmi_2018_sum_stats.txt > bmi_tmp_maf.txt # Filter out SNPs with MAF < 0.01
	# (2b) P filter
awk '$9<0.05' bmi_tmp_maf.txt > bmi_tmp_maf_p.txt # Filter out SNPs with P >= 0.05
	# (2c) Bi-allelic filter
awk '$4=="G" || $4=="A" || $4=="C" ||$4=="T"' bmi_tmp_maf_p.txt > bmi_tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$5=="G" || $5=="A" || $5=="C" ||$5=="T"' bmi_tmp_bialleles.txt > bmi_tmp_maf_p_biallelic.txt # Filter out multi-allelic non-effect alleles
awk 'BEGIN { print "CHR\tPOS\tSNP\tTested_Allele\tOther_Allele\tFreq_Tested_Allele_in_HRS\tBETA\tSE\tP\tN" } { print }' bmi_tmp_maf_p_biallelic.txt > bmi_tmp_maf_p_biallelic_header.txt # Re-add header to file
cp bmi_tmp_maf_p_biallelic_header.txt "$TMPDIR"/bmi/output_bmi
mv "$TMPDIR"/bmi/output_bmi/bmi_tmp_maf_p_biallelic_header.txt "$TMPDIR"/bmi/output_bmi/processed_bmi_2018_sum_stats.txt
	# (2d) 1KGP filter
# TO-DO: Check against 1KGP reference SNPs
# TO-DO: Write output file > ${output_file}
	# (2e) DEPICT formatting
awk 'BEGIN { OFS = "\t" } FNR>1 {print $3,$1,$2,$9}' ${output_file} > bmi_tmp_depict.txt # Extract columns 'SNP', 'CHR', 'BP' and 'P'
awk 'BEGIN { print "SNP\tChr\tPos\tP" }{ print }' bmi_tmp_depict.txt > ${depict_file} # Rename columns to 'SNP', 'Chr', 'Pos', 'P' for final DEPICT file


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l bmi_2018_sum_stats.txt >> ${summary_file} # Original file size
wc -l bmi_tmp_maf.txt >> ${summary_file} # Effect of MAF filter
wc -l bmi_tmp_maf_p.txt >> ${summary_file} # Effect of P filter
wc -l bmi_tmp_maf_p_biallelic_header.txt >> ${summary_file} # Effect of bi-allelic filter
wc -l ${output_file} >> ${summary_file} # Processed file size
wc -l ${depict_file} >> ${summary_file} # DEPICT file size
	# (3b) Export output directory
cp -r "$TMPDIR"/bmi/output_bmi $HOME/Koen/GWAS_data
	# (3c) Clean temporary files
rm bmi_tmp*
