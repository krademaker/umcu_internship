#!/usr/bin/python

# GOAL:     Restructure scRNA-seq expression data to serve as input for DEPICT tissue enrichment analysis.
# INPUT:    COUNT DATA: raw gene expression count data (see README for details)
#           CLUSTER DATA: clustering identity of cells (see README for details)
#           MOUSE MAPPING: mouse gene ids mapped to mouse gene names (see README for details)
#           MOUSE-HUMAN MAPPING: mouse gene ids mapped to human ortholog gene ids (see README for details)
# AUTHOR:   Koen Rademaker
# DATE:     23 April 2019

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

########## Set paths to files ##########
filename_count_data = '{}/FILENAME'.format(path)
filename_cluster_data = '{}/FILENAME'.format(path)
filename_output = '{}/FILENAME'.format(path)
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

########## Average cells per cluster ##########
count_data_cluster_averaged = get_average_by_celltype(count_data_normalized, cluster_data)

########## Map mouse genes to human ##########
count_data_cluster_averaged_human, unmapped = to_ensembl(mouse_mapping, mouse_human_mapping, count_data_cluster_averaged)

########## Standardize across clusters ##########
count_data_cluster_averaged_human_standardized = standardize(count_data_cluster_averaged_human)
count_data_cluster_averaged_human_standardized.sort_index(axis=1, inplace=True)

########## Save result to output file ##########
count_data_cluster_averaged_human_standardized.to_csv(filename_output, index=True, header=True, sep='\t')
