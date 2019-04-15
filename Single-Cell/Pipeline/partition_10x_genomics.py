#!/usr/bin/python3

# ABOUT: Script to partition full 10x Genomics 1.3 million mouse brain cells dataset for easier processing
# REQUIRED: - HDF5 gene / cell matrix (see README)
# AUTHOR: Koen Rademaker
# DATE: 15 April 2019


# Step 1 - Import packages
import sys
import scanpy as sc
import pandas as pd
import numpy as np


# Step 2 - Set Scanpy settings
sc.logging.print_versions()
sc.settings.verbosity=5
sc.settings.logfile='partition_10x_genomics.log'
sc.settings.max_memory=60
sc.settings.n_jobs=4


# Step 3 - Set variables
filename_input_h5=sys.argv[1]
filename_cluster_ids=sys.argv[2]
genome='mm10'
n_cells=1306217
n_partitions=4
partition_size=int(np.round(n_cells/n_partitions, decimals=0))


# Step 4 - Load data
sc_data=sc.read_10x_h5(filename_input_h5, genome)
sc_data.var_names_make_unique()
sc_data.obs['cell_labels'] = pd.read_csv(filename_cluster_ids, header=None, skiprows=1, dtype='category')[1].values
sc.logging.print_memory_usage()


# Step 5 - Partition data and save as H5AD files
partition=sc_data[sc_data.obs_names[0:partition_size], :]
partition.write_h5ad('partition_1.h5ad')

partition=sc_data[sc_data.obs_names[partition_size:2*partition_size], :]
partition.write_h5ad('partition_2.h5ad')

partition=sc_data[sc_data.obs_names[2*partition_size:3*partition_size], :]
partition.write_h5ad('partition_3.h5ad')

partition=sc_data[sc_data.obs_names[3*partition_size:4*partition_size+1], :]
partition.write_h5ad('partition_4.h5ad')
