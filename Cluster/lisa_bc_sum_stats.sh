#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process BC (Breast Cancer summary statistics data with (1) MAF > 0.01, (2) P < 0.05, (3) bi-allelic variants, (4) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Michailidou et al. (2017) BC summary statistics, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 27 May 2019


# (1) File organization
work_dir="$TMPDIR"/ea # Create work directory
mkdir ${work_dir}
out_dir="$TMPDIR"/ea/output # Create output directory
mkdir ${out_dir}

log_file=${out_dir}/log_bc.txt
output_file=${out_dir}/processed_bc_2017_sum_stats.txt
depict_file=${out_dir}/depict_bc_2017.txt

cp $HOME/Koen/GWAS_data/sum_stats/bc_2017_sum_stats.txt ${work_dir} # Copy summary statistics to work directory
cp $HOME/1K/g1000_eur.synonyms ${work_dir} # Copy SNP synonyms to work directory
cd ${work_dir} # Move to work directory


# (2) Data processing
    # (2a) Removal of SNPs without rs ID
cut -f1-10 bc_2017_sum_stats.txt > tmp_cohort # Select only specific cohort data
grep 'rs' tmp_cohort > tmp_cohort_rs # Select only SNPs with rs ID
awk -v n=1 '{s = gensub("^(([^:]*:){"n"}).*$", "\\1", 1); sub(".$","",s); print $0, s}' tmp_cohort_rs > tmp_cohort_rs_col # Reformat rs ID
sed -i 's/\r//' tmp_cohort_rs_col # Reformat file
sed -e 's/ /\t/g' tmp_cohort_rs_col > tmp
mv tmp tmp_cohort_rs_col
cut -f10 tmp_cohort_rs_col > tmp_cohort_rs_col_1 # Rearrange columns, starting with SNP, CHR and BP
cut -f2-9 tmp_cohort_rs_col > tmp_cohort_rs_col_2
paste tmp_cohort_rs_col_1 tmp_cohort_rs_col_2 > tmp_cohort_rearranged
    # (2b) MAF filter
awk '$6>0.01 && $6<0.99' tmp_cohort_rearranged > tmp_maf.txt # Filter out SNPs with MAF < 0.01
	# (2c) P filter
awk '$9<0.05' tmp_maf.txt > tmp_maf_p.txt # Filter out SNPs with P < 0.05
	# (2d) Bi-allelic filter
awk '$4=="G" || $4=="A" || $4=="C" || $4=="T"' tmp_maf_p.txt > tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$5=="G" || $5=="A" || $5=="C" || $5=="T"' tmp_bialleles.txt > tmp_maf_p_biallelic.txt # Filter out multi-allelic non-effect alleles
  # (2e) 1KGP filter
awk 'OFS=""; FNR>2 { for (i=1; i<=NF; i++) print "rs",$i }' g1000_eur.synonyms > tmp_synonyms.txt # Reformat SNP synonyms to include 'rs' prefix
mv tmp_synonyms.txt g1000_eur.synonyms
awk 'NR == FNR{c[$1]++;next};c[$1] > 0' g1000_eur.synonyms tmp_maf_p_biallelic.txt > tmp_maf_p_biallelic_1kgp.txt # Filter out SNPs not included in 1KGP list
awk 'BEGIN { OFS="\t"; print "SNP\tCHR\tBP\tA1\tA2\tMAF\tBETA\tSE\tP" } {print $0}' tmp_maf_p_biallelic_1kgp.txt > ${output_file} # Re-add header to file
	# (2f) DEPICT formatting
awk 'BEGIN { OFS="\t"; print "SNP\tChr\tPos\tP" } FNR>1 {print $1,$2,$3,$9}' ${output_file} > ${depict_file} # Extract SNP, Chr, Pos and P columns


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l bc_2017_sum_stats.txt >> ${log_file} # Original file size
wc -l tmp_maf.txt >> ${log_file} # Effect of MAF filter
wc -l tmp_maf_p.txt >> ${log_file} # Effect of P filter
wc -l tmp_maf_p_biallelic.txt >> ${log_file} # Effect of bi-allelic filter
wc -l ${output_file} >> ${log_file} # Processed file size
wc -l ${depict_file} >> ${log_file} # DEPICT file size
	# (3b) Export output directory files
cp -r ${out_dir}/* $HOME/Koen/GWAS_data/processed
	# (3c) Clean temporary files
rm tmp*
