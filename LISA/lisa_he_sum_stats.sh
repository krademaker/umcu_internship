#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process height summary statistics data with (1) MAF > 0.01, (2) P < 0.05, (3) bi-allelic variants, (4) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Yengo et al. (2018) height summary statistics, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 30 March 2019


# (1) File organization
work_dir="$TMPDIR"/he # Create work directory
mkdir ${work_dir}
out_dir="$TMPDIR"/he/output # Create output directory
mkdir ${out_dir}

log_file=${out_dir}/log_he.txt
output_file=${out_dir}/processed_he_2018_sum_stats.txt
depict_file=${out_dir}/depict_he_2018.txt

cp $HOME/Koen/GWAS_data/sum_stats/he_2018_sum_stats.txt ${work_dir} # Copy summary statistics to work directory
cp $HOME/1K/g1000_eur.synonyms ${work_dir} # Copy SNP synonyms to work directory
cd ${work_dir} # Move to work directory
head -n 1 he_2018_sum_stats.txt > tmp_header.txt


# (2) Data processing
	# (2a) MAF filter
awk '$6>0.01 && $6<0.99' he_2018_sum_stats.txt > tmp_maf.txt # Filter out SNPs with MAF < 0.01
	# (2b) P filter
awk '$9<0.05' tmp_maf.txt > tmp_maf_p.txt # Filter out SNPs with P>= 0.05
	# (2c) Bi-allelic filter
awk '$4=="G" || $4=="A" || $4=="C" ||$4=="T"' tmp_maf_p.txt > tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$5=="G" || $5=="A" || $5=="C" ||$5=="T"' tmp_bialleles.txt > tmp_maf_p_biallelic.txt # Filter out multi-allelic non-effect alleles
cat tmp_header.txt tmp_maf_p_biallelic.txt > tmp_cat.txt # Re-add header to file
mv tmp_cat.txt tmp_maf_p_biallelic.txt
  # (2d) 1KGP filter
awk 'BEGIN { OFS=""; print "SNP" } FNR>2 { for (i=1; i<=NF; i++) print "rs",$i }' g1000_eur.synonyms > tmp_synonyms.txt # Reformat SNP synonyms to include 'rs' prefix
mv tmp_synonyms.txt g1000_eur.synonyms
awk 'NR == FNR{c[$1]++;next};c[$3] > 0' g1000_eur.synonyms tmp_maf_p_biallelic.txt > ${output_file} # Filter out SNPs not included in 1KGP list
	# (2e) DEPICT formatting
awk 'BEGIN { OFS = "\t" } FNR>1 {print $3,$1,$2,$9}' ${output_file} > tmp_depict.txt # Extract columns 'SNP', 'CHR', 'BP' and 'P'
awk 'BEGIN { print "SNP\tChr\tPos\tP" }{ print }' tmp_depict.txt > ${depict_file} # Reformat header for DEPICT file


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l he_2018_sum_stats.txt >> ${log_file} # Original file size
wc -l tmp_maf.txt >> ${log_file} # Effect of MAF filter
wc -l tmp_maf_p.txt >> ${log_file} # Effect of P filter
wc -l tmp_maf_p_biallelic.txt >> ${log_file} # Effect of bi-allelic filter
wc -l ${output_file} >> ${log_file} # Processed file size
wc -l ${depict_file} >> ${log_file} # DEPICT file size
	# (3b) Export output directory files
cp -r ${out_dir}/* $HOME/Koen/GWAS_data/processed
	# (3c) Clean temporary files
rm tmp*
