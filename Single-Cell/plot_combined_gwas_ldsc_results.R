# TITLE:    plot_combined_gwas_ldsc_results.R
# ABOUT:    Script to plot LDSC results and highlight significant cell types for different GWASs.
# INPUT:    file_path_10x_genomics_lvl1: Combined LDSC results for 10x Genomics dataset level 1 (n=16)
# AUTHOR:   Koen Rademaker
# DATE:     29 July 2019


########## Load required libraries ##########
library(ggplot2)
library(dplyr)


########## Set variables ##########
gwas_order <- c('SCZ (2018)', 'EA (2018)', 'MDD (2018)', 'Intelligence (2018)', 'BMI (2015)', 'Breast cancer (2017)')
cell_type_order <- c('Astrocytes (1)', 'Astrocytes (2)', 'Oligodendrocytes', 'Microglia', 'Enteric neurons', 'Enteric glial cells', 'Endothelial cells', 'Vascular endothelial cells', 'Cajal-Retzius cells', 'Interneurons (1)', 'Interneurons (2)', 'Neurons (glutamatergic)', 'Neurons', 'Intermediate progenitors', 'Neuroblasts (1)', 'Neuroblasts (2)')
file_path_10x_genomics_lvl1 <- '~/umcu_internship/Single-Cell/data/10x_Genomics_16_cell_types_combined_GWAS_LDSC_results_heatmap.csv'


########## Plot LDSC results for 10x Genomics, k-means clustering (k=20), 16 cell types ##########
# (1) Load data
heatmap_10x_genomics <- read.csv(file_path_10x_genomics_lvl1)
# (2) Mark significant P-values (Bonferroni-threshold: 0.05/(16*6) )
significant_frames = heatmap_10x_genomics[heatmap_10x_genomics$Significant, c('GWAS', 'Cell.type')]
significant_frames <- significant_frames %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'MDD (2018)', which(gwas_order == 'MDD (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'Intelligence (2018)', which(gwas_order == 'Intelligence (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'BMI (2015)', which(gwas_order == 'BMI (2015)'))) %>% 
    mutate(Cell.type = as.character(Cell.type)) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Endothelial cells', which(cell_type_order == 'Endothelial cells'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Interneurons (1)', which(cell_type_order == 'Interneurons (1)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Interneurons (2)', which(cell_type_order == 'Interneurons (2)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Neuroblasts (1)', which(cell_type_order == 'Neuroblasts (1)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Neuroblasts (2)', which(cell_type_order == 'Neuroblasts (2)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Neurons', which(cell_type_order == 'Neurons'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Neurons (glutamatergic)', which(cell_type_order == 'Neurons (glutamatergic)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Oligodendrocytes', which(cell_type_order == 'Oligodendrocytes')))
significant_frames <- data.frame(lapply(significant_frames,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_10x_genomics) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order), fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'springgreen4', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames, size = 0.75, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(14) +
    theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
    theme(title = element_text(size = 9), axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))
