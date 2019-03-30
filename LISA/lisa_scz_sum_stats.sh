#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process schizophrenia summary statistics data with (1) INFO > 0.9, (2) MAF > 0.01, (3) P < 0.05, (4) bi-allelic variants, (5) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Pardiñas et al. (2018) summary statistics & high-quality imputed SNPs, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 30 March 2019


# (1) File organization
work_dir="$TMPDIR"/scz # Create work directory
mkdir ${work_dir}
out_dir="$TMPDIR"/scz/output # Create output directory
mkdir ${out_dir}

log_file=${out_dir}/log_scz.txt
output_file=${out_dir}/processed_scz_2018_sum_stats.txt
depict_file=${out_dir}/depict_scz_2018.txt

cp $HOME/Koen/GWAS_data/sum_stats/scz_2018_sum_stats.txt ${work_dir} # Copy summary statistics to work directory
cp $HOME/Koen/GWAS_data/sum_stats/scz_2018_hq_imputed_snps.txt ${work_dir} # Copy high quality imputed SNPs to work directory
cp $HOME/1K/g1000_eur.synonyms ${work_dir} # Copy SNP synonyms to work directory
cd ${work_dir} # Move to work directory
head -n 1 scz_2018_sum_stats.txt > tmp_header.txt


# (2) Data processing
  # (2a) MAF filter
awk '$2>0.01 && $2<0.99' scz_2018_sum_stats.txt > tmp_maf.txt # Filter out SNPs with MAF < 0.01
  # (2b) INFO filter
awk '{print $1}' scz_2018_hq_imputed_snps.txt > tmp_snps.txt # Extract column 'SNP'
mv tmp_snps.txt scz_2018_hq_imputed_snps.txt
awk 'NR == FNR{c[$1]++;next};c[$1] > 0' scz_2018_hq_imputed_snps.txt tmp_maf.txt > tmp_maf_info.txt # Merge columns 'SNP' where values overlap in both files
  # (2c) P filter
awk '$9<0.05' tmp_maf_info.txt > tmp_maf_info_p.txt # Filter out SNPs with P >= 0.05
  # (2d) Bi-allelic filter
awk '$5=="G" || $5=="A" || $5=="C" ||$5=="T"' tmp_maf_info_p.txt > tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$6=="G" || $6=="A" || $6=="C" ||$6=="T"' tmp_bialleles.txt > tmp_maf_info_p_biallelic.txt # Filter out multi-allelic non-effect alleles
cat tmp_header.txt tmp_maf_info_p_biallelic.txt > tmp_cat.txt # Re-add header to file
mv tmp_cat.txt tmp_maf_info_p_biallelic.txt
  # (2e) Beta calculation
awk '{print $0, "\tBETA"}' tmp_header.txt > tmp_beta_header.txt # Update header
mv tmp_beta_header.txt tmp_header.txt
awk 'FNR>1 {print $0, log($7)}' tmp_maf_info_p_biallelic.txt > tmp_beta.txt # Calculate beta as log(OR)
cat tmp_header.txt tmp_beta.txt > tmp_cat.txt # Re-add header to file
mv tmp_cat.txt tmp_maf_info_p_biallelic_beta.txt
  # (2f) rsID filter
grep 'rs' tmp_maf_info_p_biallelic_beta.txt > tmp_maf_info_p_biallelic_beta_rs.txt # Filter out SNPs without rsID
cat tmp_header.txt tmp_maf_info_p_biallelic_beta_rs.txt > tmp_cat.txt # Re-add header to file
mv tmp_cat.txt > tmp_maf_info_p_biallelic_beta_rs.txt
  # (2g) 1KGP filter
awk 'BEGIN { OFS=""; print "SNP" } FNR>2 { for (i=1; i<=NF; i++) print "rs",$i }' g1000_eur.synonyms > tmp_synonyms.txt # Reformat SNP synonyms to include 'rs' prefix
mv tmp_synonyms.txt g1000_eur.synonyms
awk 'NR == FNR{c[$1]++;next};c[$2] > 0' g1000_eur.synonyms tmp_maf_info_p_biallelic_beta_rs.txt > ${output_file} # Filter out SNPs not included in 1KGP list
	# (2h) DEPICT formatting
awk 'BEGIN { OFS = "\t" } FNR>1 {print $1,$3,$4,$9}' ${output_file} > tmp_depict.txt # Extract columns 'SNP', 'CHR', 'BP' and 'P'
awk 'BEGIN { print "SNP\tChr\tPos\tP" }{ print }' tmp_depict.txt > ${depict_file} # Reformat header for DEPICT file


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l scz_2018_sum_stats.txt >> ${log_file} # Original file size
wc -l tmp_maf.txt >> ${log_file} # Effect of MAF filter
wc -l tmp_maf_info.txt >> ${log_file} # Effect of INFO filter
wc -l tmp_maf_info_p.txt >> ${log_file} # Effect of P filter
wc -l tmp_maf_info_p_biallelic.txt >> ${log_file} # Effect of bi-allelic filter
wc -l tmp_maf_info_p_biallelic_beta_rs.txt >> ${log_file} # Effect of rsID filter
wc -l ${output_file} >> ${log_file} # Processed file size
wc -l ${depict_file} >> ${log_file} # DEPICT file size
	# (3b) Export output directory files
cp -r ${out_dir}/* $HOME/Koen/GWAS_data/processed
	# (3c) Clean temporary files
rm tmp*
