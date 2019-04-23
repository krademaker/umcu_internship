#!/usr/bin/python3

# ABOUT: Script for QC on partitions of 10x Genomics 1.3 million mouse brain cells data
# REQUIRED: - H5AD annotated gene / cell matrix
# ARGUMENTS:    1 - Path to H5AD annotated gene / cell matrix (see README)
#               2 - Mean + 3*SD limit for % mtDNA / cell
#               3 - Mean - 3*SD limit for UMI count / cell
#               4 - Mean + 3*SD limit for UMI count / cell
#               5 - Defined miminum for gene count / cell
#               6 - Mean + 3*SD limit for gene count / cell
# AUTHOR: Koen Rademaker
# DATE: 18 April 2019


# Step 1 - Import packages
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys


# Step 2 - Set Scanpy settings
sc.logging.print_versions()
sc.settings.verbosity = 5
sc.settings.logfile=sys.argv[1]+'_qc.log'
sc.settings.max_memory=60


# Step 3 - Set variables
filename_input_h5ad=sys.argv[1]
filename_output_matrix=sys.argv[1]+'_matrix.txt'
mtdna_max=float(sys.argv[2])
umi_count_min=10**(float(sys.argv[3]))
umi_count_max=10**(float(sys.argv[4]))
gene_count_min=int(sys.argv[5])
gene_count_max=float(sys.argv[6])


# Step 4 - Load data
sc_data=sc.read_h5ad(filename_input_h5ad)
sc.logging.msg('Running QC for '+filename_input_h5ad)
sc.logging.print_memory_usage()


# Step 5 - Run QC on the data
# Remove cells with no UMI counts
sc.logging.msg('Running QC on UMI counts / cell >= 1 '+str(sc_data.shape))
sc.pp.filter_cells(sc_data, min_counts=1)
sc.logging.print_memory_usage()
# Remove cells with a log10 UMI count 3 SDs below the mean
sc.logging.msg('Running QC on log10 UMI counts / cell > -3 SD  '+str(sc_data.shape))
sc.pp.filter_cells(sc_data, min_counts=umi_count_min)
sc.logging.print_memory_usage()
# Remove cells with a log10 UMI count 3 SDs above the mean
sc.logging.msg('Running QC on log10 UMI counts / cell < +3 SD '+str(sc_data.shape))
sc.pp.filter_cells(sc_data, max_counts=umi_count_max)
sc.logging.print_memory_usage()
# Remove cells with a gene count below the defined minimum
sc.logging.msg('Running QC on gene counts / cell > defined minimum '+str(sc_data.shape))
sc.pp.filter_cells(sc_data, min_genes=gene_count_min)
sc.logging.print_memory_usage()
# Remove cells with a gene count 3 SDs above the mean
sc.logging.print_memory_usage()
sc.logging.msg('Running QC on gene counts / cell < +3 SD '+str(sc_data.shape))
sc.pp.filter_cells(sc_data, max_genes=gene_count_max)
# Remove cells with a % mtDNA 3 SDs above the mean
sc.logging.msg('Running QC on %mtDNA / cell < +3 SD '+str(sc_data.shape))
sc_data=sc_data[sc_data.obs['percent_mito'] < mtdna_max,:]
sc.logging.print_memory_usage()


# Step 6 - Output H5AD file with filtered data
sc.logging.msg('Outputting updated H5AD file, final rows-columns shape: '+str(sc_data.shape))
sc_data.write_h5ad(filename_input_h5ad)


# Step 7 - Transpose and output data to gene / cell matrix
sc_data.T.to_df().to_csv(filename_output_matrix, index=True, header=True, sep='\t', compression='gzip')
