#!/bin/bash

# TITLE:	calculate_ld_scores.sh
# GOAL:		Calculate LD scores for specificity deciles of cell types
# INPUT:	LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types
# INPUT:	1000 Genomes Project phase 3 plink files, HapMap3 SNPs (will automatically be downloaded and organized)
# OUTPUT:	LD score output files; .M (total number of SNPs), .M_5_50 (number of SNPs with MAF > 5%), .ldscore.gz (SNP LD scores)
# AUTHOR:	Koen Rademaker
# DATE:		3 July 2019


########## Initialize script ##########
# Set folder paths
ldsc_dir=~/umcu_internship/Single-Cell/LDSC/ldsc
files_dir=~/umcu_internship/Single-Cell/LDSC/Files
annotation_dir=${files_dir}/Annotation
ld_scores_dir=${files_dir}/LD-Scores
tmp_dir=${files_dir}/tmp
# Create folder paths
mkdir -p ${tmp_dir}
mkdir -p ${ld_scores_dir}
# Declare variables
declare -a cell_types=("Glutamatergic_neurons" "Neuroblasts_1" "Astrocytes_1" "Neuroblasts_2" "Intermediate_progenitors" "Enteric_glial_cells" "Interneurons_1" "Neurons" "Interneurons_2" "Vascular_endothelial_cells" "Enteric_neurons" "Cajal_Retzius_cells" "Oligodendrocytes" "Astrocytes_2" "Endothelial_cells" "Microglia")
declare -a specificity_deciles=(N 1 2 3 4 5 6 7 8 9 10)


########## Organize LDSC ##########
# Download background data
echo "-- Downloading 1000 Genomes Project phase 3 plink files--"
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz -P ${files_dir}
tar -xvzf ${files_dir}/1000G_Phase3_plinkfiles.tgz
echo "-- Downloading HapMap3 SNPs --"
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz -P ${files_dir}
tar -xvzf ${files_dir}/hapmap3_snps.tgz
# Declare variables
phase3_1000G_dir=${files_dir}/1000G_EUR_Phase3_plink
hapmap3_dir=${files_dir}/hapmap3_snps
ld_window=1
dataset="10x_Genomics"
# Activate LDSC environment
conda activate ldsc


########## Calculate LD scores with annotation files ##########
# (1) Iterate over cell types
# (2) Iterate over chromosomes
# (3) Copy annotation file to temporary folder
# (4) Calculate LD scores
# (5) Merge files together for export
# (6) Export
# (7) Clear temporary folder
