#!/usr/bin/python


# ABOUT: Script for QC and restructuring of 10x Genomics mouse brain scRNA-seq data
# REQUIRED: - HFD5 gene / cell matrix (see README)
#           - Clustering output file (see README)
#           - ENSEMBL to MGI mapping file (see README)
#           - ENSEMBL mouse to human mapping file (see README)
# AUTHOR: Koen Rademaker
# DATE: 11 April 2019


# Step 1 - Import packages
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import os
import sys


# Step 2 - Import 'perslab-sc-library' repository
path_to_script = os.path.realpath(__file__).split('/')
path = '/'.join(path_to_script[0:(len(path_to_script)-2)])
sys.path.append('{}/Single-Cell/perslab-sc-library'.format(path))
from dropseq import get_average_by_celltype, standardize, normalize
from gene_mapping import to_ensembl


# Step 3 - Set Scanpy settings
sc.logging.print_versions()
sc.settings.verbosity = 5
sc.settings.logfile='lisa_10x_1M_qc_restructure.log'
sc.settings.autosave = True
sc.settings.max_memory=40
sc.settings.n_jobs=4
sc.settings.set_figure_params(dpi=300)


# Step 4 - Set variables
filename_input_h5=sys.argv[1]
filename_cluster_ids=sys.argv[2]
filename_mouse_mapping=sys.argv[3]
filename_mouse_to_human=sys.argv[4]
genome='mm10'
filename_h5ad_out='10x_1M_neurons_qc.h5ad'
filename_restructured_out='10x_1M_neurons_qc_restructured.csv'
filename_mapped_out='10x_1M_neurons_qc_human_mapped.csv'


# Step 5 - Load data
sc_data=sc.read_10x_h5(filename_input_h5, genome)
sc_data.var_names_make_unique()
sc_data.obs['cell_labels'] = pd.read_csv(filename_cluster_ids, header=None, skiprows=0, dtype='category')[1].values


# Step 6 - Load mouse to human gene mapping
mouse_mapping=pd.read_csv(filename_mouse_mapping, compression='gzip', index_col=1, sep='\t')
mouse_mapping=mouse_mapping[[isinstance(x, str) for x in mouse_mapping.index]]
mouse_to_human=pd.read_csv(filename_mouse_to_human, compression='gzip', index_col=0, sep='\t')


# Step 7 - Determine cell quality metrics (percentage mitochondrial genes, n.o. counts, n.o. genes)
mito_genes=sc_data.var_names.str.startswith('mt-')
sc_data.obs['percent_mito'] = np.sum(sc_data[:, mito_genes].X, axis=1).A1 / np.sum(sc_data.X, axis=1).A1
sc_data.obs['n_counts'] = sc_data.X.sum(axis=1).A1
sc.pp.filter_cells(sc_data, min_genes=0)


# Step 8 - Plot pre-QC metrics
sc.pl.violin(sc_data, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save='pre_qc_metrics.png')
sc.pl.scatter(sc_data, x='n_counts', y='percent_mito', save='pre_qc_counts_mito.png')
sc.pl.scatter(sc_data, x='n_counts', y='n_genes', save='pre_qc_counts_genes.png')


# Step 9 - Run QC on the data
sc.pp.filter_cells(sc_data, min_counts=1) # Retain cells with >1 UMI counts
sc.pp.filter_cells(sc_data, min_genes=500) # Retain cells with >500 genes
sc.pp.filter_cells(sc_data, max_genes=4000) # Retain cells with <4000 genes
sc.pp.filter_genes(sc_data, min_cells=1) # Retain genes occuring in >1 cell
sc.pp.filter_genes(sc_data, min_counts=1) # Retain genes with >1 UMI counts
sc_data.write_h5ad(filename_h5ad_out) # Save H5AD object for later use
sc_data=sc_data[sc_data.obs['percent_mito'] < 0.1,:] # Retain cells with <10% mitochondrial genes


# Step 10 - Plot post-QC metrics
sc.pl.violin(sc_data, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save='post_qc_metrics.png')
sc.pl.scatter(sc_data, x='n_counts', y='percent_mito', save='post_qc_counts_mito.png')
sc.pl.scatter(sc_data, x='n_counts', y='n_genes', save='post_qc_counts_genes.png')


# Step 11 - Restructure data
normalized_df=normalize(sc_data.to_df()) # Normalize gene expression across cells
cell_cluster_ids=pd.read_csv(filename_cluster_ids, header=None, skiprows=1, index_col=0) # Load cluster IDs for cells
cell_cluster_ids.columns=['cluster_id'] # Set DataFrame column name
averaged_per_cell_type_df=get_average_by_celltype(normalized_df.T, cell_cluster_ids) # Average gene expression across cell types
mapped_to_human,not_mapped=to_ensembl(mouse_mapping, mouse_to_human, averaged_per_cell_type_df) # Map to human genes
standardized_by_cell_type_df=standardize(mapped_to_human) # Standardize gene expression across cell types
mapped_to_human.to_csv(filename_mapped_out, sep='\t') # Export mouse-to-human mapping results
standardized_by_cell_type_df.to_csv(filename_restructured_out, sep='\t') # Export final restructured expression matrix
