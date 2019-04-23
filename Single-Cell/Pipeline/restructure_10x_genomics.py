#!/usr/bin/python


# ABOUT: Script to restructure 10x Genomics mouse brain scRNA-seq data for DEPICT tissue enrichment analysis
# REQUIRED: - Gene / cell matrix (see README)
# AUTHOR: Koen Rademaker
# DATE: 23 April 2019


# Step 1 - Import packages
import numpy as np
import os
import pandas as pd
import sys


# Step 2 - Import 'perslab-sc-library' repository
path_to_script = os.path.realpath(__file__).split('/')
path = '/'.join(path_to_script[0:(len(path_to_script)-2)])
sys.path.append('{}/scripts/perslab-sc-library'.format(path))
from dropseq import get_average_by_celltype, standardize, normalize
from gene_mapping import to_ensembl


# Step 3 - Set variables
filename_count_data = sys.argv[1]
filename_cluster_data = sys.argv[2]
filename_mouse_mapping = sys.argv[3]
filename_mouse_to_human = sys.argv[4]
filename_output='10x_genomics_post_qc_restructured_ENSEMBL_v96_GRCh37.txt.gz'


# Step 4 - Load and format input data
count_data = pd.read_csv(filename_count_data, compression='gzip', sep='\t', index_col=0)

cluster_data = pd.read_csv(filename_cluster_data, sep=',', index_col=0, skiprows=0)
cluster_data.columns = ['cluster_id']

mouse_mapping = pd.read_csv(filename_mouse_mapping, compression='gzip', index_col=1, sep='\t')
mouse_mapping = mouse_mapping[[isinstance(x, str) for x in mouse_mapping.index]]

mouse_human_mapping = pd.read_csv(filename_mouse_human_mapping, compression='gzip', index_col=0, sep='\t')


# Step 5 - Restructure data
# Normalize count data
count_data_normalized = normalize(count_data)
del count_data
# Average cells per cluster
count_data_cluster_averaged = get_average_by_celltype(count_data_normalized, cluster_data)
del count_data_normalized
# Map mouse genes to human
count_data_cluster_averaged_human, unmapped = to_ensembl(mouse_mapping, mouse_human_mapping, count_data_cluster_averaged)
del count_data_cluster_averaged
# Standardize across clusters
count_data_cluster_averaged_human_standardized = standardize(count_data_cluster_averaged_human)
count_data_cluster_averaged_human_standardized.sort_index(axis=1, inplace=True)
del count_data_cluster_averaged_human


# Step 6 - Save result to output file
count_data_cluster_averaged_human_standardized.to_csv(filename_output, index=True, header=True, sep='\t', compression='gzip')
