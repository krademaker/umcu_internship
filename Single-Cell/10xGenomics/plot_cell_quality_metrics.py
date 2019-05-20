#!/usr/bin/python

# ABOUT:        Script to plot pre- and post-QC cell quality metrics
# PARAMETERS:   - H5AD object with pre-qc data (sc_data_pre_qc)
#               - H5AD object with post_qc data (sc_data_post_qc)
#               - Number of cells, used in plot titles (n_cells)
# AUTHOR:	Koen Rademaker
# DATE:		17 May 2019


########## Import packages and settings ##########
import pandas as pd
import numpy as np
import scanpy as sc
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams


########## Set settings ##########
matplotlib.use('Agg')
rcParams['figure.figsize'] = 6, 6
rcParams['savefig.dpi'] = 600
matplotlib.style.use('seaborn-bright')


########## Load data ##########
sc_data_pre_qc = sc.read_h5ad(sys.argv[1])
sc_data_post_qc = sc.read_h5ad(sys.argv[2])
n_cells = sys.argv[3]


########## Calculate mean and standard deviation ##########
mtdna_mean = np.mean(sc_data_pre_qc.obs['percent_mito'])*100
mtdna_sd = np.std(sc_data_pre_qc.obs['percent_mito'])*100

umi_mean = np.mean(np.log10(sc_data_pre_qc.obs['n_counts']))
umi_sd = np.std(np.log10(sc_data_pre_qc.obs['n_counts']))

gene_mean = np.mean(sc_data_pre_qc.obs['n_genes'])
gene_sd = np.std(sc_data_pre_qc.obs['n_genes'])

cell_mean = np.mean(sc_data_pre_qc.var['n_cells'])
cell_sd = np.std(sc_data_pre_qc.var['n_cells'])


########## Generate pre-QC plots ##########
# Calculate percentage of cell below and above the UMI cut-off values.
perc_below_min_umi = np.round((len(np.log10(sc_data_pre_qc.obs['n_counts']).where(np.log10(sc_data_pre_qc.obs['n_counts']) < umi_mean-3*umi_sd).dropna())/len(np.log10(sc_data_pre_qc.obs['n_counts'])))*100, decimals=2)
perc_above_max_umi = np.round((len(np.log10(sc_data_pre_qc.obs['n_counts']).where(np.log10(sc_data_pre_qc.obs['n_counts']) > umi_mean+3*umi_sd).dropna())/len(np.log10(sc_data_pre_qc.obs['n_counts'])))*100, decimals=2)
# Plot annotated UMI distribution
plt.hist(np.log10(sc_data_pre_qc.obs['n_counts']), bins=100, color='blue')
plt.xlabel('UMIs per cell (log10)')
plt.ylabel('Frequency')
plt.title('UMI distribution for '+n_cells+' cells')
plt.tight_layout()
plt.axvline(x=umi_mean-3*umi_sd, color='red', linestyle='--', label='< µ-3σ (~'+str(perc_below_min_umi)+'% of cells)')
plt.axvline(x=umi_mean, color='black', linestyle='--', label='µ')
plt.axvline(x=umi_mean+3*umi_sd, color='red', linestyle=':', label='> µ+3σ (~'+str(perc_above_max_umi)+'% of cells)')
plt.legend(loc='upper right', fontsize='small')
plt.show()
plt.savefig('PRE_QC_umi_distribution_annotated.png', dpi=600)
plt.close()

# Calculate percentage of cell below and above the gene cut-off values.
perc_below_min_genes=np.round((len(sc_data_pre_qc.obs['n_genes'].where(sc_data_pre_qc.obs['n_genes'] < 1000).dropna())/len(sc_data_pre_qc.obs['n_genes']))*100, decimals=2)
perc_above_max_genes=np.round((len(sc_data_pre_qc.obs['n_genes'].where(sc_data_pre_qc.obs['n_genes'] > gene_mean+3*gene_sd).dropna())/len(sc_data_pre_qc.obs['n_genes']))*100, decimals=2)
# Plot annotated gene distribution
plt.hist(sc_data_pre_qc.obs['n_genes'], bins=100, color='blue')
plt.xlabel('Genes per cell')
plt.ylabel('Frequency')
plt.title('Gene distribution for '+n_cells+' cells')
plt.tight_layout()
plt.axvline(x=1000, color='red', linestyle='--', label='< 1000 (~'+str(perc_below_min_genes)+'% of cells)')
plt.axvline(x=gene_mean, color='black', linestyle='--', label='µ')
plt.axvline(x=gene_mean+3*gene_sd, color='red', linestyle=':', label='> µ+3σ (~'+str(perc_above_max_genes)+'% of cells)')
plt.legend(loc='upper right', fontsize='small')
plt.show()
plt.savefig('PRE_QC_gene_distribution_annotated.png', dpi=600)
plt.close()

# Plot annotated %mtDNA distribution
sc_data_pre_qc.obs['percent_mito'] = sc_data_pre_qc.obs['percent_mito'].apply(lambda x: x*100)
perc_exceeding_mito_threshold=np.round((len(sc_data_pre_qc.obs['percent_mito'].where(sc_data_pre_qc.obs['percent_mito'] > mtdna_mean+3*mtdna_sd).dropna())/len(sc_data_pre_qc.obs['percent_mito']))*100, decimals=2)

plt.hist(sc_data_pre_qc.obs['percent_mito'], bins=100, range=[0, 20], color='blue')
plt.xlabel('% mtDNA per cell')
plt.ylabel('Frequency')
plt.title('% mtDNA distribution for '+n_cells+' cells')
plt.tight_layout()
plt.axvline(x=mtdna_mean, color='black', linestyle='--', label='µ')
plt.axvline(x=mtdna_mean+3*mtdna_sd, color='red', linestyle=':', label='> µ+3σ (~'+str(perc_exceeding_mito_threshold)+'% of cells)')
plt.legend(loc='upper right', fontsize='small')
plt.show()
plt.savefig('PRE_QC_mtdna_distribution_annotated.png', dpi=600)
plt.close()


########## Generate post-QC plots ##########
# Plot UMI distribution
plt.hist(np.log10(sc_data_post_qc.obs['n_counts']), bins=100, color='red')
plt.xlabel('UMIs per cell (log10)')
plt.ylabel('Frequency')
plt.title('UMI distribution for +'n_cells'+ cells after QC')
plt.tight_layout()
plt.show()
plt.savefig('POST_QC_umi_distribution.png', dpi=600)
plt.close()

# Plot gene distribution
plt.hist(sc_data_post_qc.obs['n_genes'], bins=100, color='red')
plt.xlabel('Genes per cell')
plt.ylabel('Frequency')
plt.title('Gene distribution for '+n_cells+' cells after QC')
plt.tight_layout()
plt.show()
plt.savefig('POST_QC_gene_distribution.png', dpi=600)
plt.close()

# Plot %mtDNA distribution
sc_data_post_qc.obs['percent_mito'] = sc_data_post_qc.obs['percent_mito'].apply(lambda x: x*100)
plt.hist(sc_data_post_qc.obs['percent_mito'], bins=100, range=[0, 20], color='red')
plt.xlabel('% mtDNA per cell')
plt.ylabel('Frequency')
plt.title('% mtDNA distribution for '+n_cells+' cells after QC')
plt.tight_layout()
plt.show()
plt.savefig('POST_QC_mtdna_distribution.png', dpi=600)
plt.close()
