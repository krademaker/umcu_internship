#!/bin/bash
#SBATCH -n 4
#SBATCH -t 48:00:00
#SBATCH -p normal
#SBATCH --mem=60G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl

# GOAL:     Pipeline for QC and restructuring of the full 10x Genomics mouse brain scRNA-seq dataset
# INPUT:    See the specific steps for the required input data
# AUTHOR:   Koen Rademaker
# DATE:     23 April 2019

########## Assign variables to script paths ##########
partition_script=$HOME/Koen/SC_data/10x_Genomics/scripts/partition_10x_genomics.py
get_qc_script=$HOME/Koen/SC_data/10x_Genomics/scripts/get_qc_metrics_10x_genomics.py
merge_qc_script=$HOME/Koen/SC_data/10x_Genomics/scripts/merge_qc_metrics_10x_genomics.py
qc_script=$HOME/Koen/SC_data/10x_Genomics/scripts/qc_10x_genomics.py
restructure_script=$HOME/Koen/SC_data/10x_Genomics/scripts/restructure_10x_genomics.py

########## Set variables for file and directory paths ##########
work_dir="$TMPDIR"/10x_genomics
out_dir="$TMPDIR"/10x_genomics/output
h5_file="$TMPDIR"/10x_genomics/1M_neurons_filtered_gene_bc_matrices_h5.h5
cluster_file="$TMPDIR"/10x_genomics/clusters.csv
mouse_mapping_file="$TMPDIR"/10x_genomics/ensembl_v96_ensembl_genename_Mm.txt.gz
mouse_human_mapping_file="$TMPDIR"/10x_genomics/ensembl_v96_Mm_Hs_GRCh37.txt.gz

########## Create directories ##########
mkdir ${work_dir}
mkdir ${out_dir}

########## Copy files to work directory ##########
cp $HOME/Koen/SC_data/10x_Genomics/1M_neurons_filtered_gene_bc_matrices_h5.h5 ${work_dir}
cp $HOME/Koen/SC_data/10x_Genomics/analysis/clustering/graphclust/clusters.csv ${work_dir}
cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v96_ensembl_genename_Mm.txt.gz ${work_dir}
cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v96_Mm_Hs_GRCh37.txt.gz ${work_dir}

########## Move to work directory ##########
cd ${work_dir}

########## Execute Python script to partition gene expression data into partitions (n=4) ##########
# ARGUMENTS:    ${h5_file}: 10x Genomics 1.3 million mouse brain cells data as cell / gene (1306127 X 27999) matrix in H5-format.
#               ${cluster_file}: Clustering identities for all cells, derived from graph-based clustering, n=60.
# OUTPUT:       partition_1.h5ad: Partition of 326532 cells from original cell / gene matrix, including clustering annotation, in H5AD format.
#               partition_2.h5ad: Partition of 326532 cells from original cell / gene matrix, including clustering annotation, in H5AD format.
#               partition_3.h5ad: Partition of 326532 cells from original cell / gene matrix, including clustering annotation, in H5AD format.
#               partition_4.h5ad: Partition of 326531 cells from original cell / gene matrix, including clustering annotation, in H5AD format.
$HOME/Koen/Python-3.6.0/python ${partition_script} ${h5_file} ${cluster_file}


########## Execute Python script to get QC metrics (% mtDNA, UMI counts, gene counts) for each partition ##########
# ARGUMENTS:    ${partition_file}: One of four aforementioned H5AD files.
# OUTPUT:       ${partition_file}_mtDNA_distribution.png: Plot with the distribution of % mtDNA / cell
#               ${partition_file}_UMI_distribution.png: Plot with the distribution of log10 UMI counts / cell
#               ${partition_file}_Gene_distribution.png: Plot with the distribution of gene counts / cell
#               ${partition_file}_percent_mito.csv: Values for % mtDNA / cell, in CSV format.
#               ${partition_file}_n_counts.csv: Values for log10 UMI counts / cell, in CSV format.
#               ${partition_file}_n_genes.csv: Values for gene counts / cell, in CSV format.
for partition_file in partition_1.h5ad partition_2.h5ad partition_3.h5ad partition_4.h5ad
do
    $HOME/Koen/Python-3.6.0/python ${get_qc_script} ${partition_file}
done

########## Set variables for output files from ${get_qc_script} ##########
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

########## Execute Python script to merge QC metrics for all partitions ##########
# ARGUMENTS:    ${mtdna_out_1}: Aforementioned % mtDNA / cell file for partition 1.
#               ${umi_out_1}: Aforementioned log10 UMI counts / cell file for partition 1.
#               ${gene_out_1}: Aforementioned gene counts / cell file for partition 1.
#               ${mtdna_out_2}: Aforementioned % mtDNA / cell file for partition 2.
#               ${umi_out_2}: Aforementioned log10 UMI counts / cell file for partition 2.
#               ${gene_out_2}: Aforementioned gene counts / cell file for partition 2.
#               ${mtdna_out_3}: Aforementioned % mtDNA / cell file for partition 3.
#               ${umi_out_3}: Aforementioned log10 UMI counts / cell file for partition 3.
#               ${gene_out_3}: Aforementioned gene counts / cell file for partition 3.
#               ${mtdna_out_4}: Aforementioned % mtDNA / cell file for partition 4.
#               ${umi_out_4}: Aforementioned log10 UMI counts / cell file for partition 4.
#               ${gene_out_4}: Aforementioned gene counts / cell file for partition 4.
# OUTPUT:       merged_qc_params_10x_genomics.log: File containing the mean and standard deviation for the QC metrics
$HOME/Koen/Python-3.6.0/python ${merge_qc_script} ${mtdna_out_1} ${mtdna_out_2} ${mtdna_out_3} ${mtdna_out_4} ${umi_out_1} ${umi_out_2} ${umi_out_3} ${umi_out_4} ${gene_out_1} ${gene_out_2} ${gene_out_3} ${gene_out_4}

########## Execute Python script to run QC for each partition ##########
# ARGUMENTS:    ${partition_file}: One of four aforementioned H5AD files.
#               Max. % mtDNA / cell limit, set to ~ 8.6%.
#               Min. log10 UMI counts / cell limit, set to ~ 2.95.
#               Max. log10 UMI counts / cell limit, set to ~ 4.30.
#               Min. gene count / cell limit, set to 1000.
#               Max. gene count / cell limit, set to ~ 4335.
# OUTPUT:       Updated ${partition_file} with low-quality cells removed.
#               ${partition_file}_matrix.txt: Gene / cell matrix of post-QC results for partition, in TSV format.
for partition_file in partition_1.h5ad partition_2.h5ad partition_3.h5ad partition_4.h5ad
do
    $HOME/Koen/Python-3.6.0/python ${qc_script} ${partition_file} 0.08658564886 2.964166266 4.295002666 1000 4335.257054968
done

########## Merging post-QC gene / cell matrices into a single matrix ##########
post_qc_matrix="$TMPDIR"/10x_genomics/post_qc_matrix.txt
gunzip partition_*.h5ad_matrix.txt.gz
paste -d'\t' partition_1.h5ad_matrix.txt <(cut -f2- partition_2.h5ad_matrix.txt) <(cut -f2- partition_3.h5ad_matrix.txt) <(cut -f2- partition_4.h5ad_matrix.txt) | column -s'\t' > ${post_qc_matrix}
gzip partition_*.h5ad_matrix.txt.gz
gzip ${post_qc_matrix}
post_qc_matrix="$TMPDIR"/10x_genomics/post_qc_matrix.txt.gz

########## Execute Python script to restructure full gene / cell matrix ##########
# ARGUMENTS:    ${post_qc_matrix}: Aforementioned post-QC gene / cell matrix.
#               ${cluster_file}: Aforementioned clustering identities file.
#               ${mouse_mapping_file}: ENSEMBL v96 mapping of mouse gene names to mouse gene IDs (GRCm38).
#               ${mouse_human_mapping_file}: ENSEMBL v96 mapping of mouse genes (GRCm38) to human orthologs (GRCh37).
# OUTPUT:       ${restructured_matrix}: Restructured output in the form of cluster-level relative mean expression gene / cell matrix,
$HOME/Koen/Python-3.6.0/python ${restructure_script} ${post_qc_matrix} ${cluster_file} ${mouse_mapping_file} ${mouse_human_mapping_file}
restructured_matrix="$TMPDIR"/10x_genomics_post_qc_restructured_ENSEMBL_v96_GRCh37.txt.gz

########## Export QC metric output files (logs, text, images) from scratch to HOME ##########
mv *.log ${out_dir}
mv *.png ${out_dir}
mv partition_*.csv ${out_dir}
cp ${out_dir}/merged_qc_params_10x_genomics.log $HOME/Koen/SC_data/10x_Genomics/output/qc_metrics
cp ${out_dir}/*.csv $HOME/Koen/SC_data/10x_Genomics/output/qc_metrics
cp ${out_dir}/*.png $HOME/Koen/SC_data/10x_Genomics/output/qc_metrics
cp ${out_dir}/*.log $HOME/Koen/SC_data/10x_Genomics/output/logs

########## Export H5AD files from scratch to HOME ##########
mv *.h5ad ${out_dir}
cp ${out_dir}/*.h5ad $HOME/Koen/SC_data/10x_Genomics/output/partitions

########## Export matrices (partitions, post-QC, restructured) from scratch to HOME ##########
mv partition_*.h5ad_matrix.txt.gz ${out_dir}
mv ${post_qc_matrix} ${out_dir}
mv ${restructured_matrix} ${out_dir}
cp ${out_dir}/*.txt.gz $HOME/Koen/SC_data/10x_Genomics/output/matrices

########## Deprecated, to be removed in the future but kept for reference ##########
# cp $HOME/Koen/SC_data/10x_Genomics/output/partitions/partition* ${work_dir} # Copy partitioned H5AD files to work directory
# cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v83_ensembl_mgisymbol_Mm.tab.gz ${work_dir} # Copy mouse mapping to work directory
# mouse_mapping_file="$TMPDIR"/10x_genomics/ensembl_v83_ensembl_mgisymbol_Mm.tab.gz
# cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v82_Mm_Hs.tab.gz ${work_dir} # Copy mouse-to-human mapping to work directory
# mouse_human_mapping_file="$TMPDIR"/10x_genomics/ensembl_v82_Mm_Hs.tab.gz
