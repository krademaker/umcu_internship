#!/bin/bash
  
# TITLE:        partition_h2.sh
# GOAL:         Partition heritability (h2) to specificity deciles of cell types.
# INPUT:        LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types, LDSC LD score files 
# INPUT:	1000 Genomes Project phase 3 baseline model LD scores, regressio weights, allele frequencies (see README)
# OUTPUT:       Partitioned h2 output files
# AUTHOR:       Koen Rademaker
# DATE:         5 July 2019


########## Initialize script ##########
# Set paths
ldsc_dir=~/Git/umcu_internship/Single-Cell/LDSC/ldsc
files_dir=~/Git/umcu_internship/Single-Cell/LDSC/Files
sum_stats_dir=${files_dir}/sum_stats
weights_dir=${files_dir}/weights_hm3_no_hla
ref_ld_dir=${files_dir}/ref_ld
baseline_dir=${files_dir}/1000G_EUR_Phase3_baseline
frq_dir=${files_dir}/1000G_Phase3_frq
out_dir=${files_dir}/partitioned_h2
# Declare variables
declare -a cell_types=("Microglia" "Endothelial_cells" "Astrocytes_2" "Oligodendrocytes" "Cajal_Retzius_cells" "Enteric_neurons" "Vascular_endothelial_cells" "Interneurons_2" "Neurons" "Interneurons_1" "Enteric_glial_cells" "Intermediate_progenitors" "Neuroblasts_2" "Astrocytes_1" "Neuroblasts_1" "Glutamatergic_neurons")


########## Organize LDSC ##########
# Download background data
echo "-- Downloading regression weights --"
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz -P ${files_dir}
tar -xvzf ${files_dir}/weights_hm3_no_hla.tgz
echo "-- Downloading baseline model LD scores --"
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz -P ${files_dir}
tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
echo "-- Downloading allele frequencies --"
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz -P ${files_dir}
tar -xvzf 1000G_Phase3_frq.tgz
# Declare variables
sum_stats_file="SCZ.sumstats.gz"
# Activate LDSC environment


########## Partition heritability across cell type specificity deciles ##########
python ldsc.py 
	--h2 ${sum_stats_dir}/${sum_stats_file} \
	--w-ld-chr ${weights_dir}/weights. \
	--ref-ld-chr ${baseline_dir}/baseline. \	# ALL MODELS
	--overlap-annot \
	--frqfile-chr ${frq_dir}/1000G.EUR.QC. \
	--out ${out_dir}/SCZ \
	--print-coefficients
