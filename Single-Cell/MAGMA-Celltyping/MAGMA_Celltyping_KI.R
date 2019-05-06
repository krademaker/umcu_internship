# TITLE:    MAGMA_Celltyping_KI.R
# ABOUT:    Script to run MAGMA Celltyping analysis on Karolinska Institue (KI) mouse brain scRNA-seq data
# INPUT:    ctd_allKI.rda: RData object containing specificity and mean expression data
#           SCZ GWAS
#           BMI GWAS
#           MDD GWAS
#           HEIGHT GWAS
#           EA GWAS
#           IQ GWAS
# AUTHOR:   Koen Rademaker
# DATE:     1 May 2019

########## Load required libraries ##########
library(One2One)
library(MAGMA.Celltyping)

########## Set parameters ##########
genome_ref_dir = '~/Git/umcu_internship/Single-Cell/MAGMA-Celltyping/g1000_eur'
genome_ref_path = sprintf('%s/g1000_eur', genome_ref_dir)

########## Prepare quantile groups for cell types ##########
ctd = prepare.quantile.groups(ctd_allKI,
                              specificity_species = 'mouse',
                              numberOfBins = 40)

# TODO: Later add paths for ALL summary statistic files

#################### Analysis for Pardiñas et al. 2018 schizophrenia GWAS ####################
########## Set parameters ##########
scz_gwas_sum_stats_path = ''

########## Map SNPs to genes ##########
scz_genes_out_path = map.snps.to.genes(scz_gwas_sum_stats_path,
                                       genome_ref_path = genome_ref_path)

##########

########## Linear cell type association analysis ##########
scz_ct_assoc_linear = calculate_celltype_associations(ctd,
                                                      scz_gwas_sum_stats_path,
                                                      genome_ref_path = genome_ref_path,
                                                      specificity_species = 'mouse')
scz_figs_linear = plot_celltype_associations(scz_ct_assoc_linear,
                                             ctd = ctd)

########## Top 10% cell type association analysis ##########
scz_ct_assoc_top = calculate_celltype_associations(ctd,
                                                   scz_gwas_sum_stats_path,
                                                   genome_ref_path = genome_ref_path,
                                                   EnrichmentMode = 'Top 10%')
scz_figs_top = plot_celltype_associations(scz_ct_assoc_top,
                                          ctd = ctd)

########## Merge linear results with top 10% results ##########
scz_ct_assoc_merged = merge_magma_results(scz_ct_assoc_linear, scz_ct_assoc_top)
scz_figs_merged = plot_celltype_associations(scz_ct_assoc_merged,
                                             ctd = ctd)

########## Conditional cell type association analysis ##########
scz_ct_assoc_condit = calculate_conditional_celltype_associations(ctd,
                                                                  scz_gwas_sum_stats_path,
                                                                  genome_ref_path = genome_ref_path,
                                                                  analysis_name = 'Conditional')
scz_figs_condit = plot_celltype_associations(scz_ct_assoc_condit,
                                             ctd = ctd)
