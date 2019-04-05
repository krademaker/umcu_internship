# ABOUT: Script to visualize a Manhattan- & QQ plot for GWAS summary statistics
# REQUIRED: qqman library, summary statistics fidatale
install.packages('qqman')
# AUTHOR: Koen Rademaker
# DATE: 5 April 2019

# Load qqman library
library(qqman)

# Set path to summary statistics file
sum_stats_path='./Git/umcu_internship/PLINK/sex_autosomal_sum_stats.txt'

# Load summary statistics data
sum_stats <- read.delim(sum_stats_path)
# Omit NA values from summary statistics
sum_stats <- na.omit(sum_stats)

# Create Manhattan plot, annotate rsID to SNPs with a P below 1e-5
png('gwas_sex_autosomal_manhattan_plot.png')
manhattan(sum_stats, chr='CHR', bp='BP', snp='SNP', p='P', annotatePval = 0.05)
dev.off()

# Create QQ plot
png('gwas_sex_autosomal_qq_plot.png')
qq(sum_stats$P)
dev.off()
