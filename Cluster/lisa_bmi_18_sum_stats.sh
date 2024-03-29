#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process 2018 BMI (Body Mass Index) summary statistics data with (1) MAF > 0.01, (2) P < 0.05, (3) bi-allelic variants, (4) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Yengo et al. (2018) BMI summary statistics, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATA: 27 May 2019


# (1) File organization
work_dir="$TMPDIR"/bmi # Create work directory
mkdir ${work_dir}
out_dir="$TMPDIR"/bmi/output # Create output directory
mkdir ${out_dir}

log_file=${out_dir}/log_bmi.txt
output_file=${out_dir}/processed_bmi_2018_sum_stats.txt
depict_file=${out_dir}/depict_bmi_2018.txt

cp $HOME/Koen/GWAS_data/sum_stats/bmi_2018_sum_stats.txt ${work_dir} # Copy summary statistics to work directory
cp $HOME/1K/g1000_eur.synonyms ${work_dir} # Copy SNP synonyms to work directory
cd ${work_dir} # Move to work directory
head -n 1 bmi_2018_sum_stats.txt > tmp_header.txt


# (2) Data processing
	# (2a) MAF filter
awk '$6>0.01 && $6<0.99' bmi_2018_sum_stats.txt > tmp_maf.txt # Filter out SNPs with MAF < 0.01
	# (2b) P filter
awk '$9<0.05 && $9>0' tmp_maf.txt > tmp_maf_p.txt # Filter out SNPs with P >= 0.05
	# (2c) Bi-allelic filter
awk '$4=="G" || $4=="A" || $4=="C" || $4=="T"' tmp_maf_p.txt > tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$5=="G" || $5=="A" || $5=="C" || $5=="T"' tmp_bialleles.txt > tmp_maf_p_biallelic.txt # Filter out multi-allelic non-effect alleles
  # (2d) 1KGP filter
awk 'OFS=""; FNR>2 { for (i=1; i<=NF; i++) print "rs",$i }' g1000_eur.synonyms > tmp_synonyms.txt # Reformat SNP synonyms to include 'rs' prefix
mv tmp_synonyms.txt g1000_eur.synonyms
awk 'NR == FNR{c[$1]++;next};c[$3] > 0' g1000_eur.synonyms tmp_maf_p_biallelic.txt > ${output_file} # Filter out SNPs not included in 1KGP list
cat tmp_header.txt ${output_file} > tmp_cat.txt # Re-add header to file
mv tmp_cat.txt ${output_file}
	# (2e) DEPICT formatting
awk 'BEGIN { OFS="\t"; print "SNP\tChr\tPos\tP" } FNR>1 {print $3,$1,$2,$9}' ${output_file} > ${depict_file} # Extract SNP, Chr, Pos and P columns


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l bmi_2018_sum_stats.txt >> ${log_file} # Original file size
wc -l tmp_maf.txt >> ${log_file} # Effect of MAF filter
wc -l tmp_maf_p.txt >> ${log_file} # Effect of P filter
wc -l tmp_maf_p_biallelic.txt >> ${log_file} # Effect of bi-allelic filter
wc -l ${output_file} >> ${log_file} # Processed file size
wc -l ${depict_file} >> ${log_file} # DEPICT file size
	# (3b) Export output directory files
cp -r ${out_dir}/* $HOME/Koen/GWAS_data/processed
	# (3c) Clean temporary files
rm tmp*
