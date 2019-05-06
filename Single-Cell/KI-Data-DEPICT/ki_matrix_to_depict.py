#!/usr/bin/python

# GOAL:     Restructure Karolinska Institute (KI) mouse brain scRNA-seq expression data to serve as input for DEPICT tissue enrichment analysis.
# INPUT:    COUNT DATA: cluster-level mean gene expression data (see README for details)
#           MOUSE MAPPING: mouse gene ids mapped to mouse gene names (see README for details)
#           MOUSE-HUMAN MAPPING: mouse gene ids mapped to human ortholog gene ids (see README for details)
# AUTHOR:   Koen Rademaker
# DATE:     5 May 2019

########## Import modules ##########
import numpy as np
import os
import pandas as pd
import sys

########## Load perslab code library and specific module functions ##########
path_to_script = os.path.realpath(__file__).split('/')
path = '/'.join(path_to_script[0:(len(path_to_script)-1)])
sys.path.append('{}/../perslab-sc-library'.format(path))
from gene_mapping import to_ensembl
from dropseq import standardize

########## Custom functions ##########
def save_mean_cluster_level_gene_expression(df, path):
    '''Calculate mean cluster-level gene expression and save to output file.

    :param df: DataFrame with cluster-level average gene expression, (rows: genes, columns: clusters).
    :param path: Path to write output file to.
    '''
    df['mean'] = df.iloc[:, :].mean(axis=1)
    df.to_csv(path, index=True, header=True, sep='\t')
    df.drop('mean', axis=1, inplace=True)

########## Set paths to files ##########
filename_cluster_averaged_data = sys.argv[1]
filename_output_standardized = sys.argv[2]
filename_mouse_mapping = '{}/../data/ensembl_v83_ensembl_mgisymbol_Mm.tab.gz'.format(path)
filename_mouse_human_mapping = '{}/../data/ensembl_v82_Mm_Hs.tab.gz'.format(path)
#	filename_mouse_mapping = '{}/../data/ensembl_v96_ensembl_genename_Mm.txt.gz'.format(path)
#	filename_mouse_human_mapping = '{}/../data/ensembl_v96_Mm_Hs_GRCh37.txt.gz'.format(path)
filename_output_mean_cluster_gene_expression_mouse = '{}/cluster_averaged_mouse_mean_ensembl_v83.txt'.format(path)
filename_output_mean_cluster_gene_expression_human = '{}/cluster_averaged_human_mean_ensembl_v82.txt'.format(path)


########## Load and format input data ##########
count_data_cluster_averaged = pd.read_csv(filename_cluster_averaged_data, sep='\t', index_col=0)

mouse_mapping = pd.read_csv(filename_mouse_mapping, compression='gzip', index_col=1, sep='\t')
mouse_mapping = mouse_mapping[[isinstance(x, str) for x in mouse_mapping.index]]

mouse_human_mapping = pd.read_csv(filename_mouse_human_mapping, compression='gzip', index_col=0, sep='\t')

########## Calculate and save cluster-level mean gene expression ##########
save_mean_cluster_level_gene_expression(count_data_cluster_averaged, filename_output_mean_cluster_gene_expression_mouse)

########## Map mouse genes to human, calculate and save cluster-level mean gene expression ##########
count_data_cluster_averaged_human, unmapped = to_ensembl(mouse_mapping, mouse_human_mapping, count_data_cluster_averaged)
save_mean_cluster_level_gene_expression(count_data_cluster_averaged_human, filename_output_mean_cluster_gene_expression_human)

########## Standardize across clusters ##########
count_data_cluster_averaged_human_standardized = standardize(count_data_cluster_averaged_human)

########## Save result to output file ##########
count_data_cluster_averaged_human_standardized.to_csv(filename_output_standardized, index=True, header=True, sep='\t')