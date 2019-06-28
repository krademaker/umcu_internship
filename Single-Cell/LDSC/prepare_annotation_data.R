# TITLE:	preprare_annotation_data.R
# ABOUT:	TO-DO
# INPUT:	TO-DO
# AUTHOR:	Koen Rademaker
# DATE:		28 June 2019

########## Load required libraries ##########
library(MAGMA.Celltyping)
library(One2One)
library(ggplot2)
library(dplyr)
library(tibble)

########## Declare custom functions ##########
get_specificity_quantile_genes <- function(ctd, cell_type, quantile){
    genes <- as.data.frame(ctd[[1]]$specificity_quantiles[, cell_type]) %>%
        rownames_to_column('gene') %>%
        filter(ctd_10X[[1]]$specificity_quantiles[, ct] == quantile) %>%
        pull(gene)
    return(genes)
}


########## Set variables ##########
mean_expr_path <- '~/Git/umcu_internship/Single-Cell/data/10x_Genomics_16_cell_types_mean_expression.tsv'
formal_cell_type_names <- c('Glutamatergic neurons', 'Neuroblasts 1', 'Astrocytes 1', 'Neuroblasts 2', 'Intermediate progenitors', 'Enteric glial cells', 'Interneurons 1', 'Neurons', 'Interneurons 2', 'Vascular endothelial cells', 'Enteric neurons', 'Cajal-Retzius cells', 'Oligodendrocytes', 'Astrocytes 2', 'Endothelial cells', 'Microglia')

########## Load cell type-level mean expression data & calculate specificity metric S(g,c) ##########
mean_expr_data <- read.delim(mean_expr_path, row.names = 1)
specificity_data <- mean_expr_data/rowSums(mean_expr_data)

########## Restructure ctd object for MAGMA cell type association analysis ##########
ctd_10X <- MAGMA.Celltyping::ctd_allKI
ctd_10X[[1]]$mean_exp <- mean_expr_data
ctd_10X[[1]]$specificity <- specificity_data

########## Calculate specificity quantiles in 41 bins (least-to-most specific expression) ##########
ctd_10X = prepare.quantile.groups(ctd_10X, specificity_species = 'mouse', numberOfBins = 10)
ctd_10X[[2]] <- NULL

########## Store unmapped genes ##########
unmapped_genes <- row.names(specificity_data)[!(row.names(specificity_data) %in% row.names(ctd_10X[[1]]$specificity_quantiles))]

########## Get genes for all categories ##########
for (ct in names(ctd_10X[[1]]$specificity)){
    anno_non_expr <- get_specificity_quantile_genes(ctd_10X, ct, 0)
        # get human chromosome number
        # split genes per chromosome
        # convert from gene symbol to ENSEMBL gene ID
        # save to file as 10x_Genomics.[celltype].N.[chr].Geneset
            
    
    # range 1-10
        # get human chromosome number
        # split genes per chromosome
        # convert from gene symbol to ENSEMBL gene ID
        # save to file as 10x_Genomics.[celltype].[n].[chr].Geneset
}