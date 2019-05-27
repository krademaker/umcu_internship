#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process 2015 BMI (Body Mass Index) summary statistics data with (1) MAF > 0.01, (2) P < 0.05, (3) bi-allelic variants, (4) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Locke et al. (2015) BMI summary statistics, 1KGP reference SNPs and BIM file
# AUTHOR: Koen Rademaker
# DATA: 27 May 2019


# (1) File organization
work_dir="$TMPDIR"/bmi # Create work directory
mkdir ${work_dir}
out_dir="$TMPDIR"/bmi/output # Create output directory
mkdir ${out_dir}

log_file=${out_dir}/log_bmi.txt
output_file=${out_dir}/processed_bmi_2015_sum_stats.txt
depict_file=${out_dir}/depict_bmi_2015.txt

cp $HOME/Koen/GWAS_data/sum_stats/bmi_2015_sum_stats.txt ${work_dir} # Copy summary statistics to work directory
cp $HOME/1K/g1000_eur.synonyms ${work_dir} # Copy SNP synonyms to work directory
cp $HOME/1K/g1000_eur.bim ${work_dir} # Copy BIM file to work directory
cd ${work_dir} # Move to work directory
head -n 1 bmi_2015_sum_stats.txt > tmp_header.txt


# (2) Data processing
    # (2a) Add chromosome and base-pair information
tail -n+2 bmi_2015_sum_stats.txt | cut -f 1 > tmp_rs # Extract column 'SNP'
awk 'NR == FNR{c[$1]++;next};c[$2] > 0' tmp_rs g1000_eur.bim > tmp_rs_in_g1000 # Extract SNPs occurring in both files
cut -f2,1,4 tmp_rs_in_g1000 > tmp_rs_in_g1000_chr_rs_bp # Extract columns 'CHR', 'SNP' and 'BP'
sort -k2 tmp_rs_in_g1000_chr_rs_bp > tmp_rs_in_g1000_chr_rs_bp_sorted # Sort by rs ID
sort -k1 bmi_2015_sum_stats.txt > bmi_2015_sum_stats_sorted.txt # Sort by rs ID
awk -F'\t' 'NR==FNR{c[$2]++;next};c[$1] > 0' tmp_rs_in_g1000_chr_rs_bp_sorted bmi_2015_sum_stats_sorted.txt > bmi_2015_sum_stats_sorted_in_g1000.txt # Extract SNPs occurring in both files
sort -k1 -V bmi_2015_sum_stats_sorted_in_g1000.txt > tmp # Sort summary statistics by rs ID
mv tmp bmi_2015_sum_stats_sorted_in_g1000.txt
sort -k2 -V tmp_rs_in_g1000_chr_rs_bp_sorted > tmp # Sort 'CHR', 'SNP' and 'BP' data by rs ID
mv tmp tmp_rs_in_g1000_chr_rs_bp_sorted
join -t $'\t' tmp_rs_in_g1000_chr_rs_bp_sorted bmi_2015_sum_stats_sorted_in_g1000.txt -1 2 -2 1 > bmi_2015_sum_stats_chr_bp.txt # Join two files, adding chromosome and base-pair values
	# (2b) MAF filter
awk '$6>0.01 && $6<0.99' bmi_2015_sum_stats_chr_bp.txt > tmp_maf.txt # Filter out SNPs with MAF < 0.01
	# (2c) P filter
awk '$9<0.05 && $9>0' tmp_maf.txt > tmp_maf_p.txt # Filter out SNPs with P >= 0.05
	# (2d) Bi-allelic filter
awk '$4=="G" || $4=="A" || $4=="C" || $4=="T"' tmp_maf_p.txt > tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$5=="G" || $5=="A" || $5=="C" || $5=="T"' tmp_bialleles.txt > tmp_maf_p_biallelic.txt # Filter out multi-allelic non-effect alleles
    # (2e) 1KGP filter
awk 'OFS=""; FNR>2 { for (i=1; i<=NF; i++) print "rs",$i }' g1000_eur.synonyms > tmp_synonyms.txt # Reformat SNP synonyms to include 'rs' prefix
mv tmp_synonyms.txt g1000_eur.synonyms
awk 'NR == FNR{c[$1]++;next};c[$3] > 0' g1000_eur.synonyms tmp_maf_p_biallelic.txt > tmp_maf_p_biallelic_1kgp.txt # Filter out SNPs not included in 1KGP list
awk 'BEGIN { OFS="\t"; print "SNP\tCHR\tPOS\tA1\tA2\tMAF\tB\tSE\tP\tN" } FNR>1 {print $0}' tmp_maf_p_biallelic_1kgp.txt > ${output_file} # Re-add header to file
	# (2f) DEPICT formatting
awk 'BEGIN { OFS="\t"; print "SNP\tChr\tPos\tP" } FNR>1 {print $1,$2,$3,$9}' ${output_file} > ${depict_file} # Extract SNP, Chr, Pos and P columns


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l bmi_2015_sum_stats.txt >> ${log_file} # Original file size
wc -l tmp_maf.txt >> ${log_file} # Effect of MAF filter
wc -l tmp_maf_p.txt >> ${log_file} # Effect of P filter
wc -l tmp_maf_p_biallelic.txt >> ${log_file} # Effect of bi-allelic filter
wc -l ${output_file} >> ${log_file} # Processed file size
wc -l ${depict_file} >> ${log_file} # DEPICT file size
	# (3b) Export output directory files
cp -r ${out_dir}/* $HOME/Koen/GWAS_data/processed
	# (3c) Clean temporary files
rm tmp*
