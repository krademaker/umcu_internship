#!/usr/bin/python3


# ABOUT: Script for QC on partitions of 10x Genomics 1.3 million mouse brain cells data
# REQUIRED: - H5AD annotated gene / cell matrix
# ARGUMENTS:    1 - Path to H5AD annotated gene / cell matrix
#               2 - Mean + 3*SD limit for % mtDNA / cell
#               3 - Mean - 3*SD limit for UMI count / cell
#               4 - Mean + 3*SD limit for UMI count / cell
#               5 - Mean - 3*SD limit for gene count / cell
#               6 - Mean + 3*SD limit for gene count / cell
# AUTHOR: Koen Rademaker
# DATE: 17 April 2019


# Step 1 - Import packages
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys


# Step 2 - Import 'perslab-sc-library' repository
path_to_script = os.path.realpath(__file__).split('/')
path = '/'.join(path_to_script[0:(len(path_to_script)-2)])
sys.path.append('{}/scripts/perslab-sc-library'.format(path))
from dropseq import get_average_by_celltype, standardize, normalize
from gene_mapping import to_ensembl


# Step 3 - Set Scanpy settings
sc.logging.print_versions()
sc.settings.verbosity = 5
sc.settings.logfile=sys.argv[1]+'_qc.log'
sc.settings.max_memory=60


# Step 4 - Set variables
filename_input_h5ad=sys.argv[1]
mtdna_max=float(sys.argv[2])
umi_count_min=10**(float(sys.argv[3]))
umi_count_max=10**(float(sys.argv[4]))
gene_count_min=int(sys.argv[5])
gene_count_max=float(sys.argv[6])
filename_output_h5ad='qc_'+sys.argv[1]


# Step 5 - Load data
sc_data=sc.read_h5ad(filename_input_h5ad)
sc.logging.print_memory_usage()


# Step 6 - Run QC on the data
sc.pp.filter_cells(sc_data, min_counts=1) # Remove cells with no UMI counts
sc.logging.print_memory_usage()
sc.pp.filter_cells(sc_data, min_counts=umi_count_min) # Remove cells with a log10 UMI count 3 SDs below the mean
sc.logging.print_memory_usage()
sc.pp.filter_cells(sc_data, max_counts=umi_count_max) # Remove cells with a log10 UMI count 3 SDs above the mean
sc.logging.print_memory_usage()
sc.pp.filter_cells(sc_data, min_genes=gene_count_min) # Remove cells with a gene count below the defined minimum
sc.logging.print_memory_usage()
sc.pp.filter_cells(sc_data, max_genes=gene_count_max) # Remove cells with a gene count 3 SDs above the mean
sc.logging.print_memory_usage()
sc_data=sc_data[sc_data.obs['percent_mito'] < mtdna_max,:] # Remove cells with a % mtDNA 3 SDs above the mean
sc.logging.print_memory_usage()


# Step 7 - Output H5AD file with filtered data
sc_data.write_h5ad(filename_output_h5ad)
