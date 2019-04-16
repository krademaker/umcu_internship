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
# DATE: 16 April 2019


# (0) File and variable organization
work_dir="$TMPDIR"/10x_genomics # Create work directory
mkdir ${work_dir}
out_dir="$TMPDIR"/10x_genomics/output # Create output directory
mkdir ${out_dir}

partition_script_path=$HOME/Koen/SC_data/10x_Genomics/scripts/partition_10x_genomics.py # Script to partition gene / cell matrix
get_qc_script_path=$HOME/Koen/SC_data/10x_Genomics/scripts/get_qc_metrics_10x_genomics.py # Script to get QC metrics for each partition
merge_qc_script_path=$HOME/Koen/SC_data/10x_Genomics/scripts/merge_qc_metrics_10x_genomics.py # Script to merge QC metrics for all partitions
qc_script_path=$HOME/Koen/SC_data/10x_Genomics/scripts/qc_10x_genomics.py # Script to run QC on partitions
restructure_script_path=$HOME/Koen/SC_data/10x_Genomics/scripts/restructure_10x_genomics.py # Script to restructure data

cp $HOME/Koen/SC_data/10x_Genomics/1M_neurons_filtered_gene_bc_matrices_h5.h5 ${work_dir} # Copy gene / cell matrix to work directory
h5_file="$TMPDIR"/10x_genomics/1M_neurons_filtered_gene_bc_matrices_h5.h5
cp $HOME/Koen/SC_data/10x_Genomics/analysis/clustering/graphclust/clusters.csv ${work_dir} # Copy clustering output to work directory
cluster_file="$TMPDIR"/10x_genomics/clusters.csv
#cp $HOME/Koen/SC_data/10x_Genomics/output/partition* ${work_dir} # Copy partitioned H5AD files to work directory
#   cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v83_ensembl_mgisymbol_Mm.tab.gz ${work_dir} # Copy mouse mapping to work directory
#   mouse_mapping_file="$TMPDIR"/10x_genomics/ensembl_v83_ensembl_mgisymbol_Mm.tab.gz
#   cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v82_Mm_Hs.tab.gz ${work_dir} # Copy mouse-to-human mapping to work directory
#   mouse_human_mapping_file="$TMPDIR"/10x_genomics/ensembl_v82_Mm_Hs.tab.gz
cd ${work_dir} # Move to work directory


# (1) Execute Python script for file partitioning with the following arguments:
# - Gene / cell matrix - ${h5_file}
# - Clustering output - ${cluster_file}
$HOME/Koen/Python-3.6.0/python ${partition_script_path} ${h5_file} ${cluster_file}


# (2) Execute Python script to get QC metrics for each partition with the following arguments:
# - Annotated gene / cell matrix partition - ${partition_file}
for partition_file in partition_1.h5ad partition_2.h5ad partition_3.h5ad partition_4.h5ad
do
    $HOME/Koen/Python-3.6.0/python ${get_qc_script_path} ${partition_file}
done


# (3) Update variables for new output files
mtdna_out_1=partition_1.h5ad_percent_mito.csv
mtdna_out_2=partition_2.h5ad_percent_mito.csv
mtdna_out_3=partition_3.h5ad_percent_mito.csv
mtdna_out_4=partition_4.h5ad_percent_mito.csv

umi_out_1=partition_1.h5ad_n_counts.csv
umi_out_2=partition_2.h5ad_n_counts.csv
umi_out_3=partition_3.h5ad_n_counts.csv
umi_out_4=partition_4.h5ad_n_counts.csv

gene_out_1=partition_1.h5ad_n_genes.csv
gene_out_2=partition_2.h5ad_n_genes.csv
gene_out_3=partition_3.h5ad_n_genes.csv
gene_out_4=partition_4.h5ad_n_genes.csv


# (4) Execute Python script for merging QC metrics on all partitions with the following arguments:
# - % mtDNA / cell output - ${mtdna_out_X} (n=4, one for each partition)
# - UMI counts / cell output - ${umi_out_X} (n=4, one for each partition)
# - Gene counts / cell output - ${gene_out_X} (n=4, one for each partition)
$HOME/Koen/Python-3.6.0/python ${merge_qc_script_path} ${mtdna_out_1} ${mtdna_out_2} ${mtdna_out_3} ${mtdna_out_4} ${umi_out_1} ${umi_out_2} ${umi_out_3} ${umi_out_4} ${gene_out_1} ${gene_out_2} ${gene_out_3} ${gene_out_4}


# (5) Execute Python script to get QC metrics for each partition with the following arguments:
# - Annotated gene / cell matrix partition - ${partition_file}
#for partition_file in partition_1.h5ad partition_2.h5ad partition_3.h5ad partition_4.h5ad
#do
    #$HOME/Koen/Python-3.6.0/python ${qc_script_path} ${partition_file}
#done


# TODO (6) Execute Python script for data restructuring with the following parameters:
#   $HOME/Koen/Python-3.6.0/python ${restructure_script_path} ${MERGED_FILE} ${mouse_mapping_file} ${mouse_human_mapping_file}


# (6) File export
mv *.h5ad ${out_dir} # Copy partitioned H5AD files to output folder
mv *.log ${out_dir} # Copy log files to output folder
mv *.png ${out_dir} # Copy images to output folder
mv *.csv ${out_dir} # Copy annotation files to output folder
cp -r ${out_dir}/* $HOME/Koen/SC_data/10x_Genomics/output # Copy files from scratch folder to HOME
