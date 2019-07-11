# TITLE:    plot_SCZ_combined_enrichment_results.R
# ABOUT:    Script to plot DEPICT, MAGMA and LDSC enrichment results for schizophrenia in different scRNA-seq datasets.
# INPUT:    scz_enrichment_10x: Combined schizophrenia enrichment results for the 10x Genomics dataset
# INPUT:    scz_enrichment_ki: Combined schizophrenia enrichment results for the KI dataset
# AUTHOR:   Koen Rademaker
# DATE:     24 June 2019


########## Load required libraries ##########
library(cowplot)


########## Set variables ##########
file_path_10x <- '~/umcu_internship/Single-Cell/data/10x_Genomics_SCZ_enrichment_results.csv'
file_path_ki <- '~/umcu_internship/Single-Cell/data/KI_SCZ_enrichment_results.csv'
cell_type_order_10x <- c('Neurons (glutamatergic)', 'Neuroblasts (2)', 'Neurons', 'Interneurons (1)', 'Oligodendrocytes', 'Interneurons (2)', 'Neuroblasts (1)', 'Intermediate progenitors', 'Cajal-Retzius cells', 'Vascular endothelial cells', 'Astrocytes (2)', 'Enteric glial cells', 'Enteric neurons', 'Astrocytes (1)', 'Endothelial cells', 'Microglia')
cell_type_order_ki <- c('Pyramidal (CA1)', 'Pyramidal (SS)', 'Medium spiny neuron', 'Interneurons', 'Striatal interneuron', 'Dopaminergic neuron', 'Hypoth. GABAergic neuron', 'Hypoth. DA neurons', 'Oligodendrocytes', 'Serotonergic neuron', 'Oligodendrocyte precursors', 'Embr. GABA neuron', 'DA neuroblast', 'Embr. DA neuron', 'Neuroblasts', 'Oxytocin/vasopressin neurons', 'Astrocytes / Ependymal', 'Embr. midbrain nucl. neuron', 'Radial glia like', 'Hypoth. glutamat. neuron', 'Neural progenitors', 'Vasc. leptomeningeal cells', 'Microglia', 'Endothelial-Mural')
n_gwas <- 6
threshold_10x <- -log10(0.05/(16*n_gwas))
ki_threshold <- -log10(0.05/(24*n_gwas))
ki_strict_threshold <- -log10(0.05/((24+149)*8))



########## Load schizophrenia enrichment results ##########
scz_enrichment_10x <- read.csv(file_path_10x)
scz_enrichment_ki <- read.csv(file_path_ki)


########## Calculate negative log10 of P-values ##########
scz_enrichment_10x$neg_log10_p <- -log10(scz_enrichment_10x$P)
scz_enrichment_ki$neg_log10_p <- -log10(scz_enrichment_ki$P)


########## Plot schizophrenia enrichment results for different statistical methods ##########
# (1) 10x Genomics dataset
ggplot(scz_enrichment_10x,
       aes(x = factor(Cell.type, level = cell_type_order_10x),
           y = neg_log10_p,
           fill = Method)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    geom_hline(yintercept = threshold_10x, colour = 'black') +
    theme_classic(14) +
    theme(legend.position = c(0.86, 0.9), legend.title = element_blank(), legend.text = element_text(size = 11)) +
    coord_flip() +
    ylab(expression('-log'[10] *' (P value)')) +
    xlab('')
# (2) Karolinska Institute (KI) dataset
ggplot(scz_enrichment_ki,
       aes(x = factor(Cell.type, level = cell_type_order_ki),
           y = neg_log10_p,
           fill = Method)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    geom_hline(yintercept = ki_threshold, colour = 'black') +
    geom_hline(yintercept = ki_strict_threshold, colour = 'red') +
    theme_classic(12) +
    theme(legend.position = c(0.9, 0.9), legend.title = element_blank(), legend.text = element_text(size = 11)) +
    coord_flip() +
    ylab(expression('-log'[10] *' (P value)')) +
    xlab('')
