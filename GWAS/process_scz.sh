#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process schizophrenia summary statistics data with (1) INFO > 0.9, (2) MAF > 0.01, (3) P < 0.05, (4) bi-allelic variants, (5) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Pardiñas et al. (2018) summary statistics, list of high-quality imputed (INFO > 0.9) SNPs, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 25 March 2019


# (1) File organization
summary_file="$TMPDIR"/scz/output_scz/summary.txt
output_file="$TMPDIR"/scz/output_scz/processed_scz_2018_sum_stats.txt
depict_file="$TMPDIR"/scz/output_scz/depict_scz_2018.txt
mkdir "$TMPDIR"/scz # Create work directory
mkdir "$TMPDIR"/scz/output_scz # Create output directory
cp $HOME/Koen/GWAS_data/scz_gwas/scz_2018_sum_stats.txt "$TMPDIR"/scz # Copy summary statistics file to work directory
cp $HOME/Koen/GWAS_data/scz_gwas/scz_2018_hq_imputed_snps.txt "$TMPDIR"/scz # Copy high-quality imputed SNPs file to work directory
	# TO-DO: 1KGP reference SNPs # Copy 1KGP reference SNPs file to work directory
cd "$TMPDIR"/scz # Move to work directory


# (2) Data processing
	# (2a) INFO filter
awk '{print $1}' scz_2018_hq_imputed_snps.txt > scz_tmp_snps.txt # Extract column 'SNP'
awk 'NR == FNR{c[$1]++;next};c[$1] > 0' scz_tmp_snps.txt scz_2018_sum_stats.txt > scz_tmp_info.txt # Merge columns 'SNP' where values overlap in both files
	# (2b) MAF filter`
awk '$2>0.01 && $2<0.99' scz_tmp_info.txt > scz_tmp_info_maf.txt # Filter out SNPs with MAF < 0.01
	# (2c) Bi-allelic filter
awk '$5=="G" || $5=="A" || $5=="C" ||$5=="T"' scz_tmp_info_maf.txt > scz_tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$6=="G" || $6=="A" || $6=="C" ||$6=="T"' scz_tmp_bialleles.txt > scz_tmp_info_maf_biallelic.txt # Filter out multi-allelic non-effect alleles
awk 'BEGIN{print "SNP\tFreq.A1\tCHR\tBP\tA1\tA2\tOR\tSE\tP"}{print}' scz_tmp_info_maf_biallelic.txt > scz_tmp_info_maf_biallelic_header.txt # Re-add header to file
	# (2d) Beta calculation
awk 'BEGIN { print "BETA" } FNR>1 {a=log($7); printf"%0.4f\n", a}' scz_tmp_info_maf_biallelic_header.txt > scz_tmp_beta.txt # Calculate beta as the log of column 'OR'
paste -d'\t' scz_tmp_info_maf_biallelic_header.txt scz_tmp_beta.txt > scz_tmp_info_maf_biallelic_header_beta.txt # Add column 'BETA'
	# (2e) P filter
awk '$9<0.05' scz_tmp_info_maf_biallelic_header_beta.txt > scz_tmp_info_maf_biallelic_header_beta_p.txt # Filter out SNPs with P >= 0.05
	# (2f) 1KGP filter
# TO-DO: Check against 1KGP reference SNPs
# TO-DO: Write output file > ${output_file}
	# (2g) DEPICT formatting
awk 'BEGIN { OFS = "\t" } FNR>1 {print $1,$3,$4,$9}' ${output_file} > scz_tmp_depict.txt # Extract columns 'SNP', 'CHR', 'BP' and 'P'
awk 'BEGIN { print "SNP\tChr\tPos\tP" }{ print }' scz_tmp_depict.txt > ${depict_file} # Rename columns to 'SNP', 'Chr', 'Pos', 'P' for final DEPICT file


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l scz_2018_sum_stats.txt >> ${summary_file} # Original file size
wc -l scz_2018_hq_imputed_snps.txt >> ${summary_file} # Original high-quality SNPs file size
wc -l scz_tmp_info.txt >> ${summary_file} # Effect of INFO filter
wc -l scz_tmp_info_maf.txt >> ${summary_file} # Effect of MAF filter
wc -l scz_tmp_info_maf_biallelic_header.txt >> ${summary_file} # Effect of bi-allelic filter
wc -l scz_tmp_info_maf_biallelic_header_beta_p.txt >> ${summary_file} # Effect of P filter
wc -l ${output_file} >> ${summary_file} # Processed file size
wc -l ${depict_file} >> ${summary_file} # DEPICT file size
	# (3b) Export output directory
cp -r "$TMPDIR"/scz/output_scz $HOME/Koen/GWAS_data
	# (3c) Clean temporary files
rm scz_tmp*
