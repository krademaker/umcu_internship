#!/usr/bin/python3


# ABOUT: Script to plot metrics of 10x Genomics mouse brain scRNA-seq dataset
# REQUIRED: - HFD5 gene / cell matrix (see README)
# AUTHOR: Koen Rademaker
# DATE: 15 April 2019


# Step 1 - Import packages
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')


# Step 2 - Set Scanpy settings
sc.logging.print_versions()
sc.settings.verbosity = 5
sc.settings.logfile='plot_metrics_10x_genomics.log'
sc.settings.autosave = True
sc.settings.set_figure_params(dpi=300)


# Step 3 - Set variables
filename_input_h5=sys.argv[1]
reference_genome='mm10'
perc_mito_threshold=10


# Step 4 - Import data
sc_data = sc.read_10x_h5(filename_input_h5, reference_genome)
sc_data.var_names_make_unique()


# Step 5 - Calculate metrics
mito_genes=sc_data.var_names.str.startswith('mt-')
sc_data.obs['percent_mito'] = (np.sum(sc_data[:, mito_genes].X, axis=1).A1 / np.sum(sc_data.X, axis=1).A1)*100
sc_data.obs['n_counts'] = sc_data.X.sum(axis=1).A1
sc.pp.filter_cells(sc_data, min_genes=0)
umi_mean=np.mean(np.log10(sc_data.obs['n_counts']))
umi_sd=np.std(np.log10(sc_data.obs['n_counts']))
gene_mean=np.mean(sc_data.obs['n_genes'])
gene_sd=np.std(sc_data.obs['n_genes'])
perc_exceeding_mito_threshold=np.round((len(sc_data.obs['percent_mito'].where(sc_data.obs['percent_mito'] > perc_mito_threshold).dropna())/len(sc_data.obs['percent_mito']))*100, decimals=1)


# Step 6 - Plot Scanpy metrics
sc.pl.scatter(sc_data, x='n_counts', y='n_genes', save='_counts_genes.png') # Relation between UMI counts and gene counts
sc.pl.scatter(sc_data, x='n_counts', y='percent_mito', save='_counts_mito.png') # Relation between UMI counts and % mtDNA
sc.pl.scatter(sc_data, x='n_genes', y='percent_mito', save='_genes_mito.png') # Relation between gene counts and % mtDNA


# Step 7 - Update Maplotlib settings
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['figure.figsize'] = 5, 5
rcParams['savefig.dpi'] = 300
matplotlib.style.use('classic')


# Step 8 - Plot Matplotlib metrics
plt.hist(np.log10(sc_data.obs['n_counts']), bins=20, color='blue') # Distribution of UMI counts
plt.xlabel('UMIS per cell (log10)')
plt.ylabel('Frequency')
plt.title('UMI Distribution')
plt.axvline(x=umi_mean-3*umi_sd, color='black', linestyle='--', label='-3 σ')
plt.axvline(x=umi_mean, color='red', linestyle='--', label='μ')
plt.axvline(x=umi_mean+3*umi_sd, color='black', linestyle='--', label='+3 σ')
plt.show()
plt.legend(loc='upper right')
plt.savefig('UMI_distribution.png')
plt.close()

plt.hist(sc_data.obs['n_genes'], bins=40, color='blue') # Distribution of gene counts
plt.xlabel('Genes per cell')
plt.ylabel('Frequency')
plt.title('Gene Distribution')
plt.axvline(x=gene_mean-3*gene_sd, color='black', linestyle='--', label='-3 σ')
plt.axvline(x=gene_mean, color='red', linestyle='--', label='μ')
plt.axvline(x=gene_mean+3*gene_sd, color='black', linestyle='--', label='+3 σ')
plt.legend(loc='upper right')
plt.show()
plt.savefig('Gene_distribution.png')
plt.close()

plt.hist(sc_data.obs['percent_mito'],  bins=100, range=(0, 15), color='blue') # Distribution of % mtDNA
plt.xlabel('% mtDNA per cell')
plt.ylabel('Frequency')
plt.title('% Mitochondrial DNA Distribution')
plt.axvline(x=perc_mito_threshold, color='red', linestyle='--', label='>'+str(perc_mito_threshold)+'% mtDNA ('+str(perc_exceeding_mito_threshold)+'%)')
plt.legend(loc='upper right')
plt.show()
plt.savefig('mtDNA_distribution.png')
plt.close()
