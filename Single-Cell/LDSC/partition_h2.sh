#!/bin/bash
  
# TITLE:        partition_h2.sh
# GOAL:         Partition heritability (h2) to specificity deciles of cell types.
# INPUT:        LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types, LDSC LD score files 
# INPUT:	1000 Genomes Project phase 3 baseline model LD scores, regression weights, allele frequencies (see README)
# OUTPUT:       Partitioned h2 output files
# AUTHOR:       Koen Rademaker
# DATE:         31 July 2019


########## Initialize script ##########
# Set paths
ldsc_dir=~/umcu_internship/Single-Cell/LDSC/ldsc
files_dir=~/umcu_internship/Single-Cell/LDSC/Files
sum_stats_dir=${}/Summary-Statistics
annotation_dir=${files_dir}/Annotation
ld_scores_dir=${files_dir}/LD-Scores
out_dir=${files_dir}/Partitioned-h2
miniconda_dir=~/miniconda2
# Declare variables
declare -a sum_stats=("SCZ_no_hm3" "MDD_no_hm3" "Intelligence_no_hm3" "EducationalAttainment_no_hm3" "BreastCancer_no_hm3" "BMI_no_hm3")
declare -a cell_types=("Microglia" "Endothelial_cells" "Astrocytes_2" "Oligodendrocytes" "Cajal_Retzius_cells" "Enteric_neurons" "Vascular_endothelial_cells" "Interneurons_2" "Neurons" "Interneurons_1" "Enteric_glial_cells" "Intermediate_progenitors" "Neuroblasts_2" "Astrocytes_1" "Neuroblasts_1" "Glutamatergic_neurons")
dataset="10x_Genomics"


########## Organize LDSC ##########
# Download background data
echo "-- Downloading regression weights --"
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz -P ${files_dir}
tar -xvzf ${files_dir}/1000G_Phase3_weights_hm3_no_MHC.tgz
echo "-- Downloading baseline model LD scores --"
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz -P ${files_dir}
tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
echo "-- Downloading allele frequencies --"
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz -P ${files_dir}
tar -xvzf 1000G_Phase3_frq.tgz
# Copy annotation files and LD score files to reference folder
ref_ld_dir=${files_dir}/ref_ld
mkdir -p ${ref_ld_dir}
cp ${ld_scores_dir}/*l2* ${ref_ld_dir}
cp ${annotation_dir}/*.annot.gz ${ref_ld_dir}
# Declare variables
weights_dir=${files_dir}/1000G_Phase3_weights_hm3_no_MHC
baseline_dir=${files_dir}/1000G_EUR_Phase3_baseline
frq_dir=${files_dir}/1000G_Phase3_frq
# Activate LDSC environment
source ${miniconda_dir}/bin/activate ldsc


########## Partition heritability for combined cell types across specificity deciles ##########
for gwas in ${sum_stats[@]}; do
        python ${ldsc_dir}/ldsc.py \
                --h2 ${sum_stats_dir}/${gwas}.sumstats.gz \
                --w-ld-chr ${weights_dir}/weights.hm3_noMHC. \
                --ref-ld-chr ${ref_ld_dir}/${dataset}_${cell_types[0]}.,${ref_ld_dir}/${dataset}_${cell_types[1]}.,${ref_ld_dir}/${dataset}_${cell_types[2]}.,${ref_ld_dir}/${dataset}_${cell_types[3]}.,${ref_ld_dir}/${dataset}_${cell_types[4]}.,${ref_ld_dir}/${dataset}_${cell_types[5]}.,${ref_ld_dir}/${dataset}_${cell_types[6]}.,${ref_ld_dir}/${dataset}_${cell_types[7]}.,${ref_ld_dir}/${dataset}_${cell_types[8]}.,${ref_ld_dir}/${dataset}_${cell_types[9]}.,${ref_ld_dir}/${dataset}_${cell_types[10]}.,${ref_ld_dir}/${dataset}_${cell_types[11]}.,${ref_ld_dir}/${dataset}_${cell_types[12]}.,${ref_ld_dir}/${dataset}_${cell_types[13]}.,${ref_ld_dir}/${dataset}_${cell_types[14]}.,${ref_ld_dir}/${dataset}_${cell_types[15]}.,${baseline_dir}/baseline. \
                --thin-annot \
                --overlap-annot \
                --frqfile-chr ${frq_dir}/1000G.EUR.QC. \
                --out ${out_dir}/${gwas}_combined_model \
                --print-coefficients
done

