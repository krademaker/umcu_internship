#!/usr/bin/python

# GOAL:     Restructure scRNA-seq expression data to serve as input for DEPICT tissue enrichment analysis.
# INPUT:    COUNT DATA: raw gene expression count data (see README for details)
#           CLUSTER DATA: clustering identity of cells (see README for details)
#           MOUSE MAPPING: mouse gene ids mapped to mouse gene names (see README for details)
#           MOUSE-HUMAN MAPPING: mouse gene ids mapped to human ortholog gene ids (see README for details)
# AUTHOR:   Koen Rademaker
# DATE:     3 May 2019

########## Import modules ##########
import numpy as np
import os
import pandas as pd
import sys

########## Load perslab code library and specific module functions ##########
path_to_script = os.path.realpath(__file__).split('/')
path = '/'.join(path_to_script[0:(len(path_to_script)-1)])
sys.path.append('{}/perslab-sc-library'.format(path))
from gene_mapping import to_ensembl
from dropseq import normalize, get_average_by_celltype, standardize

########## Custom functions ##########
def save_mean_cluster_level_gene_expression(df, path):
    '''Calculate mean cluster-level gene expression and save to output file.

    :param df: DataFrame with cluster-level average gene expression, (rows: genes, columns: clusters).
    :param path: Path to write output file to.
    '''
    df['mean'] = df.iloc[:, :].mean(axis=1)
    df.to_csv(path, index=True, header=True, sep='\t')

########## Set paths to files ##########
filename_count_data = '{}/Combined_plates_cleaned_Feb2019_genes1500_mit20perc.txt'.format(path)
filename_cluster_data = '{}/cluster_per_cell_allcells_20clusters.txt'.format(path)
filename_output_mean_cluster_gene_expression_mouse = '{}/JEPPE_MOUSE.txt'.format(path)
filename_output_mean_cluster_gene_expression_human = '{}/JEPPE_HUMAN.txt'.format(path)
filename_output = '{}/FILENAME.txt'.format(path)
filename_mouse_mapping = '{}/data/ensembl_v96_ensembl_genename_Mm.txt.gz'.format(path)
filename_mouse_human_mapping = '{}/data/ensembl_v96_Mm_Hs_GRCh37.txt.gz'.format(path)

########## Load and format input data ##########
count_data = pd.read_csv(filename_count_data, sep='\t', index_col=0)

cluster_data = pd.read_csv(filename_cluster_data, sep='\t', index_col=0, skiprows=0)
cluster_data.columns = ['cluster_id']

mouse_mapping = pd.read_csv(filename_mouse_mapping, compression='gzip', index_col=1, sep='\t')
mouse_mapping = mouse_mapping[[isinstance(x, str) for x in mouse_mapping.index]]

mouse_human_mapping = pd.read_csv(filename_mouse_human_mapping, compression='gzip', index_col=0, sep='\t')

########## Normalize count data (10K UMIs, log transformation, remove non-expressed genes) ##########
count_data_normalized = normalize(count_data)

########## Average cells per cluster, calculate and save cluster-level mean gene expression ##########
count_data_cluster_averaged = get_average_by_celltype(count_data_normalized, cluster_data)
save_mean_cluster_level_gene_expression(count_data_cluster_averaged, filename_output_mean_cluster_gene_expression_mouse)

########## Map mouse genes to human, calculate and save cluster-level mean gene expression ##########
count_data_cluster_averaged_human, unmapped = to_ensembl(mouse_mapping, mouse_human_mapping, count_data_cluster_averaged)
save_mean_cluster_level_gene_expression(count_data_cluster_averaged_human, filename_output_mean_cluster_gene_expression_human)

########## Standardize across clusters ##########
count_data_cluster_averaged_human_standardized = standardize(count_data_cluster_averaged_human)

########## Save result to output file ##########
count_data_cluster_averaged_human_standardized.to_csv(filename_output, index=True, header=True, sep='\t')
