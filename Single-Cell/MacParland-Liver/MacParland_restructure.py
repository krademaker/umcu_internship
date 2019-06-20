#!/usr/bin/python

# ABOUT: 	Script for MAGMA & DEPICT restructuring of MacParland et al. (2018) human liver scRNA-seq data
# REQUIRED: 	- Input gene / cell matrix in CSV format (filename_input_csv)
#		- Clustering identities of cells (filename_cluster_ids)
#		- Ensembl human gene symbol > human gene ID conversion file (filename_human_mapping)
#               - Output of MAGMA restructuring (filename_celltype_avg_out)
#		- Output of DEPICT restructuring (filename_depict_out)
# AUTHOR:	Koen Rademaker, GitHub repository 'perslab-sc-library' (https://github.com/perslab/perslab-sc-library, customized code for own purposes)
# DATE:	        20 June 2019


########## Import packages ##########
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys


########## Set variables ##########
filename_input_csv = sys.argv[1]
filename_cluster_ids = sys.argv[2]
filename_human_mapping = sys.argv[3]
filename_celltype_avg_out = sys.argv[4]
filename_depict_out = sys.argv[5]
cell_types = {1: "Hep 1", 2: "CD3+ αβ T cells", 3: "Hep 2", 4: "Inflammatory macrophages", 5: "Hep 3", 6: "Hep 4", 7: "Antibody secreting B cells", 8: "NK-like cells", 9: "γδ T cells 1", 10: "Non-inflammatory Macrophages", 11: "Periportal LSECs", 12: "Central venous LSECs", 13: "Portal endothelial cells", 14: "Hep 5", 15: "Hep 6", 16: "Mature B cells", 17: "Cholangiocytes", 18: "γδ T cells 2", 19: "Erthyroid cells", 20: "Stellate cells"}


########## Function declaration ##########
# Function to map gene symbols to gene identifiers
def to_ensembl(df_hs2hs, df):
	gene_symbols = []
	gene_ids = []
	gene_data = []
	for ix in df.index.tolist():
		if ix in df_hs2hs.index:
			if type(df_hs2hs.loc[ix, 'Ensembl Gene ID']) != str:
				for duplicate in df_hs2hs.loc[ix, 'Ensembl Gene ID']:
					gene_symbols.append(ix)
					gene_ids.append(duplicate)
					gene_data.append(df.loc[ix])
			else:
				gene_symbols.append(ix)
				gene_ids.append(df_hs2hs.loc[ix, 'Ensembl Gene ID'])
				gene_data.append(df.loc[ix])
	multi_index = [np.array(gene_symbols), np.array(gene_ids)]
	return pd.DataFrame(gene_data, index=multi_index)

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


########## Load human gene mapping data ##########
human_mapping = pd.read_csv(filename_human_mapping, compression='gzip', index_col=0, skiprows=1, names=['Ensembl Gene ID'], sep='\t')


########## Load data from CSV ##########
liver_data = sc.read_csv(filename_input_csv, first_column_names=True)
liver_data = liver_data.T
cell_labels = pd.read_csv(filename_cluster_ids, sep='\t', usecols=['CellName', 'Cluster#'], index_col=0)
cell_labels.replace({"Cluster#": cell_types}, inplace=True)
liver_data.obs['cell_labels'] = cell_labels


########## Restructure data for DEPICT ##########
cluster_id_cells = liver_data.obs['cell_labels'].to_frame()
cluster_id_cells.columns = ['cluster_id']
normalized = normalize(liver_data.to_df().T)								# Normalize cells
cell_type_averaged = get_average_by_celltype(normalized, cluster_id_cells)				# Average gene expression per cluster
cell_type_averaged.to_csv(filename_celltype_avg_out, sep='\t', header=True, index=True)                 # Export cell type-averaged data for MAGMA analysis
standardized = standardize(cluster_averaged)			 					# Standardize gene expression across cell types
standardized_ensembl = to_ensembl(human_mapping, standardized)						# Map gene symbols to Ensembl IDs
standardized_ensembl.to_csv(filename_depict_out, sep='\t', header=True, index=True)		        # Export standardized matrix for DEPICT analysis
