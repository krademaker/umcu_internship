#!/usr/bin/python3.6


# ABOUT: Script for QC / pre-processing of 10x Genomics mouse brain scRNA-seq data
# REQUIRED: 10x Genomics 1.3 million cells mouse brain dataset (20K cells subsample)
# AUTHOR: Koen Rademaker
# DATE: 8 April 2019


# Import packages
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys


# Import "perslab-sc-library" repository
path_to_script = os.path.realpath(__file__).split('/')
path = "/".join(path_to_script[0:(len(path_to_script)-2)])
sys.path.append("{}/Single-Cell/perslab-sc-library".format(path))
from dropseq import get_average_by_celltype, standardize


# Set variables
filename_input_h5='1M_neurons_neuron20k.h5'
filename_cluster_ids='20k_clusters.csv'
genome='mm10'
filename_out='1M_neurons_neuron_20k_normalized.csv'


# Load data
sc_data=sc.read_10x_h5(filename_input_h5, genome)
sc_data.var_names_make_unique()
sc_data.obs['cell_labels'] = pd.read_csv(filename_cluster_ids, header=None, dtype='category')[1].values


# Determine cell quality metrics (percentage mitochondrial genes, n.o. counts)
mito_genes=sc_data.var_names.str.startswith('mt-')
sc_data.obs['percent_mito'] = np.sum(sc_data[:, mito_genes].X, axis=1).A1 / np.sum(sc_data.X, axis=1).A1
sc_data.obs['n_counts'] = sc_data.X.sum(axis=1).A1


# Pre-process data
sc.pp.filter_cells(sc_data, min_counts=1) # Retain cells with at least 1 UMI count
sc.pp.filter_cells(sc_data, min_genes=500) # Retain cells with at least 500 genes
sc.pp.filter_genes(sc_data, min_cells=1) # Retain genes occurring in at least 1 cell
sc.pp.filter_genes(sc_data, min_counts=1) # Retain genes with at least 1 UMI count
sc_data=sc_data[sc_data.obs['percent_mito'] < 0.05,:] # Retain cells with fewer than 5% mitochondrial genes
sc_data=sc_data[sc_data.obs['n_genes'] < 4000,:] # Retain cells with fewer than 4000 genes
sc.pp.normalize_per_cell(sc_data, counts_per_cell_after=1e4, key_n_counts='original_counts') # Normalize to 10,000 counts per cell
sc.pp.log1p(sc_data) # Log transform the data


# Restructure data
sc_data_df=sc_data.to_df().T # Convert to Pandas DataFrame and transpose
cell_cluster_ids=pd.read_csv(filename_cluster_ids, sep=',', header=None, index_col=0) # Load cluster ids for cells
cell_cluster_ids.columns=['cluster_id'] # Set DataFrame column name
sc_data_avg_by_celltype=get_average_by_celltype(sc_data_df, cell_cluster_ids) # Average gene expression across cell types
sc_data_standardized_by_celltype=standardize(sc_data_avg_by_celltype) # Standardize gene expression across cell types
