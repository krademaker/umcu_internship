#!/bin/bash

#############################################
#$ -N calculate_ld_scores_10x_Genomics
#$ -S /bin/bash
#$ -cwd
#$ -e /home/hers_en/krademaker/errors/
#$ -o /home/hers_en/krademaker/output/
#$ -M K.J.Rademaker-2@umcutrecht.nl
#$ -m beas
#$ -l h_rt=1:00:00
#$ -l h_vmem=16G
#$ -pe threaded 4
#############################################

# TITLE:	hpc_calculate_ld_scores.sh
# GOAL:         Calculate LD scores for specificity deciles of cell types
# INPUT:        LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types
# INPUT:        1000 Genomes Project phase 3 plink files, HapMap3 SNPs (will automatically be downloaded and organized)
# OUTPUT:       LD score output files; .M (total number of SNPs), .M_5_50 (number of SNPs with MAF > 5%), .ldscore.gz (SNP LD scores)
# AUTHOR:       Koen Rademaker
# DATE:         3 July 2019


########## Initialize script ##########
# Set folder paths
ldsc_dir=/hpc/hers_en/krademaker/ldsc
files_dir=/hpc/hers_en/krademaker/ldsc_files
annotation_dir=${files_dir}/annotation
ld_scores_dir=${files_dir}/ld_scores
# Declare variables
declare -a cell_types=("Glutamatergic_neurons" "Neuroblasts_1" "Astrocytes_1" "Neuroblasts_2" "Intermediate_progenitors" "Enteric_glial_cells" "Interneurons_1" "Neurons" "Interneurons_2" "Vascular_endothelial_cells" "Enteric_neurons" "Cajal_Retzius_cells" "Oligodendrocytes" "Astrocytes_2" "Endothelial_cells" "Microglia")


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
source /hpc/hers_en/krademaker/miniconda2/bin/activate ldsc

########## Calculate LD scores with annotation files ##########
# (1) Iterate over cell types
for ct in ${cell_types[@]}; do
        # (2) Iterate over autosomal chromosomes
        for chr in {1..22}; do
                # (3) Calculate LD scores
                python ${ldsc_dir}/ldsc.py \
                        --l2 \
                        --bfile ${phase3_1000G_dir}/1000G.EUR.QC.${chr} \
                        --ld-wind-cm ${ld_window} \
                        --annot ${annotation_dir}/${dataset}_${ct}.${chr}.annot.gz \
                        --thin-annot \
                        --out ${ld_scores_dir}/${dataset}_${ct}.${chr} \
                        --print-snps ${hapmap3_dir}/hm.${chr}.snp
	done
done