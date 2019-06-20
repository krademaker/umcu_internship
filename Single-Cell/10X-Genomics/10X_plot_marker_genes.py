#!/usr/bin/python

# ABOUT:        Script to plot the location of marker genes in t-SNE space
# PARAMETERS:   - H5 object with gene / cell matrix (input_h5)
#               - t-SNE 2D-coordinates in comma-separated format (input_tsne)
#               - Output path to store images (out_path)
# AUTHOR:	Koen Rademaker, re-used code from http://cf.10xgenomics.com/supp/cell-exp/megacell_tutorial.html
# DATE:		27 May 2019


########## Import packages ##########
import collections
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
import tables
import sys


########## Function declaration ##########
GeneBCMatrix = collections.namedtuple('GeneBCMatrix', ['gene_ids', 'gene_names', 'barcodes', 'matrix'])

def get_matrix_from_h5(filename, genome):
    with tables.open_file(filename, 'r') as f:
        try:
            dsets = {}
            for node in f.walk_nodes('/' + genome, 'Array'):
                dsets[node.name] = node.read()
            matrix = sp_sparse.csc_matrix((dsets['data'], dsets['indices'], dsets['indptr']), shape=dsets['shape'])
            return GeneBCMatrix(dsets['genes'], dsets['gene_names'], dsets['barcodes'], matrix)
        except tables.NoSuchNodeError:
            raise Exception("Genome %s does not exist in this file." % genome)
        except KeyError:
            raise Exception("File is missing one or more required datasets.")

def get_expression(gbm, gene_name):
    gene_indices = np.where(gbm.gene_names == gene_name)[0]
    if len(gene_indices) == 0:
        raise Exception("%s was not found in list of gene names." % gene_name)
    return gbm.matrix[gene_indices[0], :].toarray().squeeze()


########## Set variables ##########
input_h5 = sys.argv[1]
input_tsne = sys.argv[2]
out_path = sys.argv[3]
matplotlib.use('Agg')
rcParams['figure.figsize'] = 6,6
rcParams['savefig.dpi'] = 1200
genome = 'mm10'
marker_genes = ['Reln', 'Ndnf', 'Stmn2', 'Tbr1', 'Gad1', 'Hes1', 'Aldoc', 'Olig1', 'Spink8', 'Igfbp7', 'Neurod6', 'Dlx2', 'Dlx6os1', 'Eomes', 'Vim', 'Cspg4', 'Pdgfra', 'Neurod1', 'Sox9', 'Id4', 'Csf1r', 'Mbp', 'Cx3cr1', 'Lefty1', 'Baiap2l1', 'Nov', 'Lpl', 'Gm26644', 'Drd1', 'Drd2', 'Ido1', 'Ppp1r1b', 'Adora2a', 'Wnt', 'Zic1', 'Chat', 'Mt3', 'Dbi', 'Zbtb20', 'Mt1', 'Mt2', 'Serpine2', 'Malat1', 'Gfap', 'Vwf', 'Olig2', 'Olig3', 'Igfbp1', 'Pcp4', 'Unc5d', 'Sox11', 'Rnd2', 'Cldn5', 'Snca', 'Zbtb20', 'H2afz', 'Hmgn2', 'Hmgb2', 'Gad1', 'Gad2', 'Nrn1', 'Pvalb', 'Mog', 'Aqp4', 'Vip', 'Snap25', 'Ptprc']


########## Load data ##########
gene_bc_matrix = get_matrix_from_h5(input_h5, genome)
tsne = pd.read_csv(input_tsne)


########## Plot marker genes ##########
for m in marker_genes:
    m_expr = get_expression(gene_bc_matrix, m.encode('utf-8'))
    plt.title(m+' expression')
    plt.scatter(tsne['TSNE-1'], tsne['TSNE-2'], c=m_expr, s=5, linewidths=0, cmap=plt.cm.Reds)
    plt.savefig(out_path+m+'_expression_heatmap.png')
    plt.close()
