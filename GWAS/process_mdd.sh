#!/bin/bash
#SBATCH -n 11
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process MDD (Major Depressive Disorder) summary statistics data with (1) INFO > 0.9, (2) MAF > 0.01, (3) P < 0.05, (4) bi-allelic variants, (5) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Wray et al. (2018) MDD summary statistics, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 25 March 2019


# (1) File organization
summary_file="$TMPDIR"/mdd/output_mdd/summary.txt
output_file="$TMPDIR"/mdd/output_mdd/processed_mdd_2018_sum_stats.txt
depict_file="$TMPDIR"/mdd/output_mdd/depict_mdd_2018.txt
mkdir "$TMPDIR"/mdd # Create work directory
mkdir "$TMPDIR"/mdd/output_mdd # Create output directory
cp $HOME/Koen/GWAS_data/control_gwas/mdd_2018_ex23andme_sum_stats.txt "$TMPDIR"/mdd # Copy summary statistics file to work directory
	# TO-DO: 1KGP reference SNPs # Copy 1KGP reference SNPs file to work directory
cd "$TMPDIR"/mdd # Move to work directory


# (2) Data processing
	# (2a) INFO filter
awk '$8>0.9' mdd_2018_ex23andme_sum_stats.txt > mdd_tmp_info.txt # Filter out SNPs with INFO <= 0.9
	# (2b) MAF filter
awk 'BEGIN { print "Freq_Avg"; OFS="\t" } FNR>1 {a=($6+$7)/2; $1=a; print $1}' mdd_tmp_info.txt > mdd_tmp_avg_maf.txt # Calculate average MAF for cases and controls
paste -d'\t' mdd_tmp_info.txt mdd_tmp_avg_maf.txt > mdd_tmp_info_avg_maf.txt # Add column 'Freq_Avg'
awk '$19>0.01 && $19<0.99' mdd_tmp_info_avg_maf.txt > mdd_tmp_info_maf.txt # Filter out SNPs with MAF < 0.01
	# (2c) Bi-allelic filter
awk '$4=="G" || $4=="A" || $4=="C" ||$4=="T"' mdd_tmp_info_maf.txt > mdd_tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$5=="G" || $5=="A" || $5=="C" ||$5=="T"' mdd_tmp_bialleles.txt > mdd_tmp_info_maf_biallelic.txt # Filter out multi-allelic non-effect alleles
awk 'BEGIN {print "CHR\tSNP\tBP\tA1\tA2\tFRQ_A_59851\tFRQ_U_113154\tINFO\tOR\tSE\tP\tngt\tDirection\tHetISqt\tHetDf\tHetPVa\tNca\tNco\tNeff\tFreq_Avg"} { print }' mdd_tmp_info_maf_biallelic.txt > mdd_tmp_info_maf_biallelic_header.txt # Re-add header to file
	# (2d) P filter
awk '$11<0.05' mdd_tmp_info_maf_biallelic_header.txt > mdd_tmp_info_maf_biallelic_header_p.txt # Filter out SNPs with P >= 0.05
	# (2e) 1KGP filter
# TO-DO: Check against 1KGP reference SNPs
# TO-DO: Write output file > ${output_file}
	# (2f) DEPICT formatting
awk 'BEGIN { OFS = "\t" } FNR>1 {print $2,$1,$3,$11}' ${output_file} > mdd_tmp_depict.txt # Extract columns 'SNP', 'CHR', 'BP' and 'P'
awk 'BEGIN { print "SNP\tChr\tPos\tP" }{ print }' mdd_tmp_depict.txt > ${depict_file} # Rename columns to 'SNP', 'Chr', 'Pos', 'P' for final DEPICT file


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l mdd_2018_ex23andme_sum_stats.txt >> ${summary_file} # Original file size
wc -l mdd_tmp_info.txt >> ${summary_file} # Effect of INFO filter
wc -l mdd_tmp_info_maf.txt >> ${summary_file} # Effect of MAF filter
wc -l mdd_tmp_info_maf_biallelic_header.txt >> ${summary_file} # Effect of bi-allelic filter
wc -l mdd_tmp_info_maf_biallelic_header_p.txt >> ${summary_file} # Effect of P filter
wc -l ${output_file} >> ${summary_file} # Processed file size
wc -l ${depict_file} >> ${summary_file} # DEPICT file size
	# (3b) Export output directory
cp -r "$TMPDIR"/mdd/output_mdd $HOME/Koen/GWAS_data
	# (3c) Clean temporary files
rm mdd_tmp*
