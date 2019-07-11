#!/bin/bash

# TITLE:	calculate_ld_scores.sh
# GOAL:		Calculate LD scores for specificity deciles of cell types
# INPUT:	LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types
# INPUT:	1000 Genomes Project phase 3 plink files, HapMap3 SNPs (will automatically be downloaded and organized)
# OUTPUT:	LD score output files; .M (total number of SNPs), .M_5_50 (number of SNPs with MAF > 5%), .ldscore.gz (SNP LD scores)
# AUTHOR:	Koen Rademaker
# DATE:		11 July 2019


########## Initialize script ##########
# Set folder paths
ldsc_dir=~/Git/umcu_internship/Single-Cell/LDSC/ldsc
files_dir=~/Git/umcu_internship/Single-Cell/LDSC/Files
annotation_dir=${files_dir}/Annotation
ld_scores_dir=${files_dir}/LD-Scores
tmp_dir=${files_dir}/tmp
# Create folder paths
mkdir -p ${tmp_dir}
mkdir -p ${ld_scores_dir}
# Declare variables
declare -a cell_types=("Microglia" "Endothelial_cells" "Astrocytes_2" "Oligodendrocytes" "Cajal_Retzius_cells" "Enteric_neurons" "Vascular_endothelial_cells" "Interneurons_2" "Neurons" "Interneurons_1" "Enteric_glial_cells" "Intermediate_progenitors" "Neuroblasts_2" "Astrocytes_1" "Neuroblasts_1" "Glutamatergic_neurons")


########## Organize LDSC ##########
# Download background data
#	echo "-- Downloading 1000 Genomes Project phase 3 plink files--"
#	wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz -P ${files_dir}
#	tar -xvzf ${files_dir}/1000G_Phase3_plinkfiles.tgz
#	echo "-- Downloading HapMap3 SNPs --"
#	wget https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz -P ${files_dir}
#	tar -xvzf ${files_dir}/hapmap3_snps.tgz
# Declare variables
phase3_1000G_dir=${files_dir}/1000G_EUR_Phase3_plink
hapmap3_dir=${files_dir}/hapmap3_snps
ld_window=1
dataset="10x_Genomics"


########## Calculate LD scores with annotation files ##########
# (1) Iterate over cell types
for ct in ${cell_types[@]}; do
	# (2) Iterate over chromosomes
	for chr in {1..22}; do
		# (3) Calculate LD scores
		python ${ldsc_dir}/ldsc.py \
			--l2 \
			--bfile ${phase3_1000G_dir}/1000G.EUR.QC.${chr} \
			--ld-wind-cm ${ld_window} \
			--annot ${annotation_dir}/${dataset}_${ct}.${chr}.annot.gz \
			--thin-annot \
			--out ${ld_scores_dir}/${dataset}_${ct}_${chr} \
			--print-snps ${hapmap3_dir}/hm.${chr}.snp
	done
done
