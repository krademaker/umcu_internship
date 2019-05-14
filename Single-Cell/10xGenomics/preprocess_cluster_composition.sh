#!/bin/bash

# GOAL: 	Compare the composition of clusters for different partitioned sizes of the 10x Genomics 1.3 million mouse brain cell dataset.
# INPUT:	N clustering files, each consisting of one cell barcode column and one cell-clustering identity column.
# AUTHOR:	Koen Rademaker
# DATE:		14 May 2019

########## Set paths to scripts and files ##########
visualization_script=~/umcu_internship/Single-Cell/10xGenomics/visualize_cluster_composition.R
cluster_file_1=/path/to/file
cluster_file_2=/path/to/file
cluster_file_N=/path/to/file

########## Calculate cluster compositions ##########
for cluster_file in ${cluster_file_1} ${cluster_file_2} ${cluster_file_N}
do
	awk -F',' '{print $2}' ${cluster_file} | sort -t ',' | uniq -c | tr -s ' ' > ${cluster_file}_composition
done

########## Visualize compositions with an R-script ##########
Rscript ${visualization_script}
