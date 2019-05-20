#!/usr/bin/python

# ABOUT: 	Script for QC, marker gene detection and DEPICT restructuring of MacParland et al. (2018) human liver scRNA-seq data
# REQUIRED: 	- Input gene / cell matrix in CSV format (filename_input_csv)
#		- Clustering identities of cells (filename_cluster_ids)
#		- Annotated gene / cell matrix with QC applied in H5AD format (filename_h5ad_qc)
#		- Output of DEPICT restructuring (filename_depict_out)
#		- Output logging file (filename_log_out)
# AUTHOR:	Koen Rademaker, GitHub repository 'perslab-sc-library' (https://github.com/perslab/perslab-sc-library, customized code for own purposes)
# DATE:		17 May 2019


########## Import packages ##########
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys


########## Set variables ##########
filename_input_csv = sys.argv[1]
filename_cluster_ids = sys.argv[2]
filename_h5ad_qc = sys.argv[3]
filename_depict_out = sys.argv[4]
filename_log_out = sys.argv[5]
cell_types = {1: "Hep 1", 2: "CD3+ αβ T cells", 3: "Hep 2", 4: "Inflammatory macrophages", 5: "Hep 3", 6: "Hep 4", 7: "Antibody secreting B cells", 8: "NK-like cells", 9: "γδ T cells 1", 10: "Non-inflammatory Macrophages", 11: "Periportal LSECs", 12: "Central venous LSECs", 13: "Portal endothelial cells", 14: "Hep 5", 15: "Hep 6", 16: "Mature B cells", 17: "Cholangiocytes", 18: "γδ T cells 2", 19: "Erthyroid cells", 20: "Stellate cells"}

########## Set Matplotlib & Scanpy settings ##########
sc.logging.print_versions()
sc.settings.verbosity = 5
sc.settings.logfile = filename_log_out


########## Function declaration ##########
# Function to normalize to 10k UMI and take log (from https://github.com/perslab/perslab-sc-library/blob/master/dropseq.py)
def normalize(df):
	dge = df.values									                        # (Koen: Replaced .as_matrix() with .values as the former is deprecated)
	col_sums = np.apply_along_axis(sum,0,dge)
	mat_dge_norm =  np.log( dge/[float(x) for x in col_sums] * 10000 + 1 )
	df_dge_norm = pd.DataFrame(mat_dge_norm,index=df.index,columns=df.columns)
	df_dge_norm.drop(df_dge_norm.index[df_dge_norm.sum(axis=1) == 0],axis=0,inplace=True)
	return df_dge_norm

# Function to average cells by cluster (from https://github.com/perslab/perslab-sc-library/blob/master/dropseq.py)
def get_average_by_celltype(df_dge,df_cluster):
	df_cluster = df_cluster.merge(df_dge.transpose(),left_index=True,right_index=True,how='inner').groupby('cluster_id',sort=False).mean().transpose()
	return df_cluster

# Function to standardize genes' expresion across cell types (from https://github.com/perslab/perslab-sc-library/blob/master/dropseq.py)
def standardize(df):
        return df.sub(df.mean(axis=1),axis=0).div(df.std(axis=1),axis=0)


########## Load data from CSV ##########
liver_data = sc.read_csv(filename_input_csv, first_column_names=True)
liver_data = liver_data.T
cell_labels = pd.read_csv(filename_cluster_ids, sep='\t', usecols=['CellName', 'Cluster#'], index_col=0)
cell_labels.replace({"Cluster#": cell_types}, inplace=True)
liver_data.obs['cell_labels'] = cell_labels


########## Calculate cell/gene metrics ##########
mito_genes=liver_data.var_names.str.startswith('MT-')
liver_data.obs['percent_mito'] = np.sum(liver_data[:, mito_genes].X, axis=1) / np.sum(liver_data.X, axis=1)
liver_data.obs['n_counts'] = liver_data.X.sum(axis=1)
sc.pp.filter_cells(liver_data, min_genes=0)
sc.pp.filter_genes(liver_data, min_cells=0)


########## Apply QC ##########
sc.pp.filter_cells(liver_data, min_counts=1500) 								# Remove cells with fewer than 1500 UMI counts
liver_data=liver_data[liver_data.obs['percent_mito'] < 0.5,:]							# Remove cells with a % mtDNA above 5%
sc.pp.filter_genes(liver_data, min_cells=3)									# Remove genes detected in less than 3 cells


########## Save data to H5AD for later use ##########
liver_data.write_h5ad(filename_h5ad_qc)


########## Restructure data for DEPICT ##########
cluster_id_cells = liver_data.obs['cell_labels'].to_frame()
cluster_id_cells.columns = ['cluster_id']
normalized = normalize(liver_data.to_df().T)									# Normalize cells
cluster_averaged = get_average_by_celltype(normalized, cluster_id_cells)					# Average gene expression per cluster
standardized = standardize(cluster_averaged)			 						# Standardize gene expression across cell types
standardized.to_csv(filename_depict_out, sep='\t', header=True, index=True)					# Export final matrix
