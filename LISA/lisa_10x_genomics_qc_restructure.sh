#!/bin/bash
#SBATCH -n 4
#SBATCH -t 6:00:00
#SBATCH -p normal
#SBATCH --mem=40G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


# ABOUT: Script for QC and restructuring of the full 10x Genomics mouse brain scRNA-seq dataset
# REQUIRED: - HFD5 gene / cell matrix (see README)
#           - Clustering output file (see README)
#           - ENSEMBL to MGI mapping file (see README)
#           - ENSEMBL mouse to human mapping file (see README)
# AUTHOR: Koen Rademaker
# DATE: 11 April 2019


# (0) File and variable organization
work_dir="$TMPDIR"/10x_genomics_qc_pp # Create work directory
mkdir ${work_dir}
out_dir="$TMPDIR"/10x_genomics_qc_pp/output # Create output directory
mkdir ${out_dir}

cp $HOME/Koen/SC_data/10x_Genomics/1M_neurons_filtered_gene_bc_matrices_h5.h5 ${work_dir} # Copy gene / cell matrix to work directory
${h5_file}="$TMPDIR"/10x_genomics_qc_pp/1M_neurons_filtered_gene_bc_matrices_h5.h5
cp $HOME/Koen/SC_data/10x_Genomics/analysis/clustering/graphclust/clusters.csv ${work_dir} # Copy clustering output to work directory
${cluster_file}="$TMPDIR"/10x_genomics_qc_pp/clusters.csv
cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v83_ensembl_mgisymbol_Mm.tab.gz ${work_dir} # Copy mouse mapping to work directory
${mouse_mapping_file}="$TMPDIR"/10x_genomics_qc_pp/ensembl_v83_ensembl_mgisymbol_Mm.tab.gz
cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v82_Mm_Hs.tab.gz ${work_dir} # Copy mouse-to-human mapping to work directory
${mouse_human_mapping_file}="$TMPDIR"/10x_genomics_qc_pp/ensembl_v82_Mm_Hs.tab.gz
cd ${work_dir} # Move to work directory


# (1) Execute Python program with the following arguments:
# - Gene / cell matrix - ${h5_file}
# - Clustering output - ${cluster_file}
# - Mouse mapping - ${mouse_mapping_file}
# - Mouse-to-human mapping - ${mouse_human_mapping_file}
script_path=$HOME/Koen/SC_data/10x_Genomics/scripts/qc_restructure_10x_genomics.py
$HOME/Koen/Python-3.6.0/python ${script_path} ${h5_file} ${cluster_file} ${mouse_mapping_file} ${mouse_human_mapping_file}


# (2) File export
mv 10x_1M_neurons_qc.h5ad ${out_dir} # Move to output directory
mv 10x_1M_neurons_qc_restructured.csv ${out_dir} # Move restructured data to output directory
mv 10x_1M_neurons_qc_human_mapped.csv ${out_dir} # Move human-mapped data to output directory
mv figures ${out_dir} # Move figures folder to output directory
cp -r ${out_dir}/* $HOME/Koen/SC_data/10x_Genomics/output # Copy files from scratch to HOME
