#!/usr/bin/python3.6


# ABOUT: Script for QC / pre-processing of 10x Genomics mouse brain scRNA-seq data
# REQUIRED: 10x Genomics 1.3 million cells mouse brain dataset (20K cells subsample)
# AUTHOR: Koen Rademaker
# DATE: 3 April 2019


# Load packages
import numpy as np
import pandas as pd
import scanpy as sc


# Set variables
H5_FILE='1M_neurons_neuron20k.h5'
GENOME='mm10'
NORMALIZED_FILE='1M_neurons_neuron_20k_normalized.csv'


# Load data
sc_data=sc.read_10x_h5(H5_FILE, GENOME)
sc_data.var_names_make_unique()


# Determine cell quality metrics (percentage mitochondrial genes, n.o. counts)
mito_genes=sc_data.var_names.str.startswith('mt-')
sc_data.obs['percent_mito'] = np.sum(sc_data[:, mito_genes].X, axis=1).A1 / np.sum(sc_data.X, axis=1).A1
sc_data.obs['n_counts'] = sc_data.X.sum(axis=1).A1


# Pre-process data
sc.pp.filter_cells(sc_data, min_counts=1) # Retain cells with at least 1 UMI count
sc.pp.filter_cells(sc_data, min_genes=200) # Retain cells with at least 200 genes
sc.pp.filter_genes(sc_data, min_cells=1) # Retain genes occurring in at least 1 cell
sc_data=sc_data[sc_data.obs['percent_mito'] < 0.05,:] # Retain cells with fewer than 5% mitochondrial genes
sc_data=sc_data[sc_data.obs['n_genes'] < 4000,:] # Retain cells with fewer than 4000 genes

sc.pp.normalize_per_cell(sc_data, counts_per_cell_after=1e4, key_n_counts='original_counts') # Normalize to 10,000 counts per cell

sc.pp.log1p(sc_data) # Log transform the data

sc_data.to_df().to_csv(NORMALIZED_FILE) # Export pre-processed, normalized gene expression matrix
