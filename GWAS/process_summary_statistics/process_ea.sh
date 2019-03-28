#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -p normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script to process EA (Educational Attainment) summary statistics data with (1) INFO > 0.9, (2) MAF > 0.01, (3) P < 0.05, (4) bi-allelic variants, (5) presence in 1000 Genomes Project (1KGP)
# REQUIRED: Lee et al. (2018) EA summary statistics, 1KGP reference SNPs
# AUTHOR: Koen Rademaker
# DATE: 29 March 2019


# (1) File organization
summary_file="$TMPDIR"/ea/output_ea/summary.txt
output_file="$TMPDIR"/ea/output_ea/processed_ea_2018_sum_stats.txt
depict_file="$TMPDIR"/ea/output_ea/depict_ea_2018.txt
mkdir "$TMPDIR"/ea # Create work directory
mkdir "$TMPDIR"/ea/output_ea # Create output directory
cp $HOME/Koen/GWAS_data/control_gwas/ea_2018_sum_stats.txt "$TMPDIR"/ea # Copy summary statistics file to work directory
cp $HOME/1K/g1000_eur.synonyms "$TMPDIR"/ea # Copy synonyms file for 1KGP SNPs to work directory
cd "$TMPDIR"/ea # Move to work directory


# (2) Data processing
	# (2a) MAF filter
awk '$6>0.01 && $6<0.99' ea_2018_sum_stats.txt > ea_tmp_maf.txt # Filter out SNPs with MAF < 0.01
	# (2b) P filter
awk '$9<0.05' ea_tmp_maf.txt > ea_tmp_maf_p.txt # Filter out SNPs with P < 0.05
	# (2c) Bi-allelic filter
awk '$4=="G" || $4=="A" || $4=="C" ||$4=="T"' ea_tmp_maf_p.txt > ea_tmp_bialleles.txt # Filter out multi-allelic effect alleles
awk '$5=="G" || $5=="A" || $5=="C" ||$5=="T"' ea_tmp_bialleles.txt > ea_tmp_maf_p_biallelic.txt # Filter out multi-allelic non-effect alleles
awk 'BEGIN { print "MarkerName\tCHR\tPOS\tA1\tA2\tEAF\tBeta\tSE\tPval" } { print }' ea_tmp_maf_p_biallelic.txt > ea_tmp_maf_p_biallelic_header.txt # Re-add header to file
  # (2d) 1KGP filter
awk 'BEGIN { OFS=""; print "SNP" } FNR>2 { for (i=1; i<=NF; i++) print "rs",$i }' g1000_eur.synonyms > ea_tmp_rs_g1000_eur.synonyms # Reformat SNP synonyms to include 'rs' prefix
awk 'NR == FNR{c[$1]++;next};c[$1] > 0' ea_tmp_rs_g1000_eur.synonyms ea_tmp_maf_p_biallelic_header.txt > ${output_file} # Filter out SNPs not included in 1KGP list
	# (2e) DEPICT formatting
awk 'BEGIN { OFS = "\t" } FNR>1 {print $1,$2,$3,$9}' ${output_file} > ea_tmp_depict.txt # Extract columns 'SNP', 'CHR', 'BP' and 'P'
awk 'BEGIN { print "SNP\tChr\tPos\tP" }{ print }' ea_tmp_depict.txt > ${depict_file} # Rename columns to 'SNP', 'Chr', 'Pos', 'P' for final DEPICT file


# (3) File export and cleaning
	# (3a) Summarize filter effect size
wc -l ea_2018_sum_stats.txt >> ${summary_file} # Original file size
wc -l ea_tmp_maf.txt >> ${summary_file} # Effect of MAF filter
wc -l ea_tmp_maf_p.txt >> ${summary_file} # Effect of P filter
wc -l ea_tmp_maf_p_biallelic_header.txt >> ${summary_file} # Effect of bi-allelic filter
wc -l ${output_file} >> ${summary_file} # Processed file size
wc -l ${depict_file} >> ${summary_file} # DEPICT file size
	# (3b) Export output directory
cp -r "$TMPDIR"/ea/output_ea $HOME/Koen/GWAS_data
	# (3c) Clean temporary files
rm ea_tmp*
