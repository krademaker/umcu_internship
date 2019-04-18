#!/usr/bin/python3

# ABOUT: Script to get QC metrics and plot them for partitions of 10x Genomics 1.3 million mouse brain cells data
# REQUIRED: - H5AD annotated gene / cell matrix
# AUTHOR: Koen Rademaker
# DATE: 18 April 2019


# Step 1 - Import packages
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams


# Step 2 - Set Scanpy settings
sc.logging.print_versions()
sc.settings.verbosity=5
sc.settings.logfile=str(sys.argv[1])+'_qc_10x_genomics.log'
sc.settings.max_memory=60


# Step 3 - Set Matplotlib settings
matplotlib.use('Agg')
rcParams['figure.figsize'] = 5, 5
rcParams['savefig.dpi'] = 300
matplotlib.style.use('classic')


# Step 4 - Set variables
filename_input_h5ad=sys.argv[1]
print(filename_input_h5ad)


# Step 5 - Load data
sc_data=sc.read_h5ad(filename_input_h5ad)


# Step 6 - Determine cell quality metrics (% mtDNA, UMI count, gene count)
# % mtDNA / cell
mito_genes=sc_data.var_names.str.startswith('mt-')
sc_data.obs['percent_mito']=np.sum(sc_data[:, mito_genes].X, axis=1).A1/np.sum(sc_data.X, axis=1).A1
sc_data.obs['percent_mito'].to_csv(str(filename_input_h5ad)+'_percent_mito.csv')
mito_mean=np.mean(sc_data.obs['percent_mito'])
mito_sd=np.std(sc_data.obs['percent_mito'])
print('Mean %mtDNA / cell: ', mito_mean)
print('SD %mtDNA / cell: ', mito_sd)
# UMI count / cell
sc_data.obs['n_counts']=sc_data.X.sum(axis=1).A1
np.log10(sc_data.obs['n_counts']).to_csv(str(filename_input_h5ad)+'_n_counts.csv')
umi_mean=np.mean(np.log10(sc_data.obs['n_counts']))
umi_sd=np.std(np.log10(sc_data.obs['n_counts']))
print('Mean UMIs / cell: ', umi_mean)
print('SD UMIs / cell: ', umi_sd)
# Gene count / cell
sc.pp.filter_cells(sc_data, min_genes=0)
sc_data.obs['n_genes'].to_csv(str(filename_input_h5ad)+'_n_genes.csv')
gene_mean=np.mean(sc_data.obs['n_genes'])
gene_sd=np.std(sc_data.obs['n_genes'])
print('Mean genes / cell: ', gene_mean)
print('SD genes / cell: ', gene_sd)


# Step 7 - Plot Matplotlib metrics
# Distribution of % mtDNA
plt.hist(sc_data.obs['percent_mito'],  bins=100, color='blue')
plt.xlabel('% mtDNA per cell')
plt.ylabel('Frequency')
plt.title('% Mitochondrial DNA Distribution')
plt.axvline(x=mito_mean-3*mito_sd, color='black', linestyle='--', label='-3 σ')
plt.axvline(x=mito_mean, color='red', linestyle='--', label='μ')
plt.axvline(x=mito_mean+3*mito_sd, color='black', linestyle='--', label='+3 σ ('+str(np.round(mito_mean+3*mito_sd,decimals=1))+')')
plt.legend(loc='upper right')
plt.show()
plt.savefig(str(sys.argv[1])+'_mtDNA_distribution.png')
plt.close()
# Distribution of UMI counts
plt.hist(np.log10(sc_data.obs['n_counts']), bins=20, color='blue')
plt.xlabel('UMIs per cell (log10)')
plt.ylabel('Frequency')
plt.title('UMI Distribution')
plt.axvline(x=umi_mean-3*umi_sd, color='black', linestyle='--', label='-3 σ')
plt.axvline(x=umi_mean, color='red', linestyle='--', label='μ')
plt.axvline(x=umi_mean+3*umi_sd, color='black', linestyle='--', label='+3 σ')
plt.show()
plt.legend(loc='upper right')
plt.savefig(str(sys.argv[1])+'_UMI_distribution.png')
plt.close()
# Distribution of gene counts
plt.hist(sc_data.obs['n_genes'], bins=40, color='blue')
plt.xlabel('Genes per cell')
plt.ylabel('Frequency')
plt.title('Gene Distribution')
plt.axvline(x=gene_mean-3*gene_sd, color='black', linestyle='--', label='-3 σ')
plt.axvline(x=gene_mean, color='red', linestyle='--', label='μ')
plt.axvline(x=gene_mean+3*gene_sd, color='black', linestyle='--', label='+3 σ')
plt.legend(loc='upper right')
plt.show()
plt.savefig(str(sys.argv[1])+'_Gene_distribution.png')
plt.close()


# Step 8 - Update H5AD file
sc_data.write_h5ad(filename_input_h5ad)
