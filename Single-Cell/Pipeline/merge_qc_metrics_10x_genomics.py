#!/usr/bin/python3

# ABOUT: Script to merge QC parameters for partitions of 10x Genomics 1.3 mouse brain cells data
# REQUIRED: - % mtDNA / cell output (n=4, one for each partition)
#           - UMI counts / cell output (n=4, one for each partition)
#           - Gene counts / cell output (n=4, one for each partition)
# AUTHOR: Koen Rademaker
# DATE: 18 April 2019


# Step 1 - Import packages
import sys
import pandas as pd
import numpy as np


# Step 2 - Set variables
filename_mtdna_1=sys.argv[1]
filename_mtdna_2=sys.argv[2]
filename_mtdna_3=sys.argv[3]
filename_mtdna_4=sys.argv[4]
filename_umi_1=sys.argv[5]
filename_umi_2=sys.argv[6]
filename_umi_3=sys.argv[7]
filename_umi_4=sys.argv[8]
filename_gene_1=sys.argv[9]
filename_gene_2=sys.argv[10]
filename_gene_3=sys.argv[11]
filename_gene_4=sys.argv[12]
filename_out='merged_qc_params_10x_genomics.log'


# Step 3 - Load data
mtdna_1=pd.read_csv(filename_mtdna_1, header=None)
mtdna_2=pd.read_csv(filename_mtdna_2, header=None)
mtdna_3=pd.read_csv(filename_mtdna_3, header=None)
mtdna_4=pd.read_csv(filename_mtdna_4, header=None)

umi_1=pd.read_csv(filename_umi_1, header=None)
umi_2=pd.read_csv(filename_umi_2, header=None)
umi_3=pd.read_csv(filename_umi_3, header=None)
umi_4=pd.read_csv(filename_umi_4, header=None)

gene_1=pd.read_csv(filename_gene_1, header=None)
gene_2=pd.read_csv(filename_gene_2, header=None)
gene_3=pd.read_csv(filename_gene_3, header=None)
gene_4=pd.read_csv(filename_gene_4, header=None)


# Step 4 - Concatenate data
# # Concatenate % mtDNA / cell data
mtdna_frames=[mtdna_1, mtdna_2, mtdna_3, mtdna_4]
mtdna_concat=pd.concat(mtdna_frames)
# Concatenate UMI count / cell data
umi_frames=[umi_1, umi_2, umi_3, umi_4]
umi_concat=pd.concat(umi_frames)
# Concatenate gene count / cell data
gene_frames=[gene_1, gene_2, gene_3, gene_4]
gene_concat=pd.concat(gene_frames)


# Step 5 - Calculate QC metrics on all cells
# Calculate mean and SD for % mtDNA / cell
mtdna_mean=np.mean(mtdna_concat[1])
mtdna_sd=np.std(mtdna_concat[1])
# Calculate mean and SD for UMI count / cell
umi_mean=np.mean(umi_concat[1])
umi_sd=np.std(umi_concat[1])
# Calculate mean and SD for gene count / cell
gene_mean=np.mean(gene_concat[1])
gene_sd=np.std(gene_concat[1])


# Step 6 - Save parameters to file
f=open(filename_out, 'w')
f.write('mtDNA mean: '+str(mtdna_mean)+'\nmtDNA SD: '+str(mtdna_sd)+'\nUMI mean: '+str(umi_mean)+'\nUMI SD: '+str(umi_sd)+'\nGene mean: '+str(gene_mean)+'\nGene SD: '+str(gene_sd))
f.close()
