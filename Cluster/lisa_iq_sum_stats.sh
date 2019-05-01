#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process IQ (Intelligence Quotient) summary statistics data with (1) INFO > 0.9, (2) MAF > 0.01, (3) P < 0.05, (4) bi-allelic variants, (5) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Savage et al. (2018) IQ summary statistics, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 30 March 2019


# (1) File organization
work_dir="$TMPDIR"/iq # Create work directory
mkdir ${work_dir}
out_dir="$TMPDIR"/iq/output # Create output directory
mkdir ${out_dir}

log_file=${out_dir}/log_iq.txt
output_file=${out_dir}/processed_iq_2018_sum_stats.txt
depict_file=${out_dir}/depict_iq_2018.txt

cp $HOME/Koen/GWAS_data/sum_stats/iq_2018_sum_stats.txt ${work_dir} # Copy summary statistics to work directory
cp $HOME/1K/g1000_eur.synonyms ${work_dir} # Copy SNP synonyms to work directory
cd ${work_dir} # Move to work directory
head -n 1 iq_2018_sum_stats.txt > tmp_header.txt


# (2) Data processing
  # (2a) MAF filter
awk '$7>0.01 && $7<0.99' iq_2018_sum_stats.txt > tmp_maf.txt # Filter out SNPs with MAF < 0.01
	# (2b) INFO filter
awk '$13>0.9' tmp_maf.txt > tmp_maf_info.txt # Filter out SNPs with INFO <= 0.9
	# (2c) P filter
awk '$11<0.05' tmp_maf_info.txt > tmp_maf_info_p.txt # Filter out SNPs with P >= 0.05
	# (2d) Bi-allelic filter
awk '$5=="g" || $5=="a" || $5=="c" || $5=="t"' tmp_maf_info_p.txt > tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$6=="g" || $6=="a" || $6=="c" || $6=="t"' tmp_bialleles.txt > tmp_maf_info_p_biallelic.txt # Filter out multi-allelic non-effect alleles
	# (2e) 1KGP filter
awk 'OFS=""; FNR>2 { for (i=1; i<=NF; i++) print "rs",$i }' g1000_eur.synonyms > tmp_synonyms.txt # Reformat SNP synonyms to include 'rs' prefix
mv tmp_synonyms.txt g1000_eur.synonyms
awk 'NR == FNR{c[$1]++;next};c[$1] > 0' g1000_eur.synonyms tmp_maf_info_p_biallelic.txt > ${output_file} # Filter out SNPs not included in 1KGP list
cat tmp_header.txt ${output_file} > tmp_cat.txt # Re-add header to file
mv tmp_cat.txt ${output_file}
	# (2f) DEPICT formatting
awk 'BEGIN { OFS="\t"; print "SNP\tChr\tPos\tP" } FNR>1 {print $1,$3,$4,$11}' ${output_file} > ${depict_file} # Extract SNP, Chr, Pos and P columns


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l iq_2018_sum_stats.txt >> ${log_file} # Original file size
wc -l tmp_maf.txt >> ${log_file} # Effect of MAF filter
wc -l tmp_maf_info.txt >> ${log_file} # Effect of INFO filter
wc -l tmp_maf_info_p.txt >> ${log_file} # Effect of P filter
wc -l tmp_maf_info_p_biallelic.txt >> ${log_file} # Effect of bi-allelic filter
wc -l ${output_file} >> ${log_file} # Processed file size
wc -l ${depict_file} >> ${log_file} # DEPICT file size
	# (3b) Export output directory files
cp -r ${out_dir}/* $HOME/Koen/GWAS_data/processed
	# (3c) Clean temporary files
rm tmp*
