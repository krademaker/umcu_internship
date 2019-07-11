# TITLE:    plot_combined_gwas_depict_results.R
# ABOUT:    Script to plot DEPICT results and highlight significant cell types for different scRNA-seq datasets.
# INPUT:    file_path_10x_genomics_lvl1: Combined GWAS DEPICT results for 10x Genomics dataset level 1 (n=16)
# INPUT:    file_path_10x_genomics_lvl2: Combined GWAS DEPICT results for 10x Genomics dataset level 2 (n=59)
# INPUT:    file_path_karolinska_lvl1: Combined GWAS DEPICT results for Karolinska Institute (KI) dataset level 1 (n=24)
# INPUT:    file_path_karolinska_lvl2: Combined GWAS DEPICT results for Karolinska Institute (KI) dataset level 2 (n=149)
# INPUT:    file_path_macparland: Combined GWAS DEPICT results for MacParland dataset
# AUTHOR:   Koen Rademaker
# DATE:     20 June 2019


########## Load required libraries ##########
library(ggplot2)
library(dplyr)


########## Set variables ##########
gwas_order <- c('SCZ (2018)', 'EA (2018)', 'MDD (2018)', 'Intelligence (2018)', 'BMI (2015)', 'Breast cancer (2017)')
cell_type_order_10x_genomics_lvl1 <- c('Astrocytes (1)', 'Astrocytes (2)', 'Oligodendrocytes', 'Microglia', 'Enteric neurons', 'Enteric glial cells', 'Endothelial cells', 'Vascular endothelial cells', 'Cajal-Retzius cells', 'Interneurons (1)', 'Interneurons (2)', 'Neurons (glutamatergic)', 'Neurons', 'Intermediate progenitors', 'Neuroblasts (1)', 'Neuroblasts (2)')
cell_type_order_10x_genomics_lvl2 <- c('Unclassified 1', 'Unclassified 2', 'Unclassified 3', 'Unclassified 4', 'Unclassified 5', 'Astrocytes 1', 'Astrocytes 2', 'Astrocytes 3', 'Astrocytes 4', 'Oligodendrocytes 1', 'Oligodendrocytes 2', 'Oligodendrocytes 3', 'Microglia', 'Radial glia 1', 'Radial glia 2', 'Radial glia 3', 'Radial glia 4', 'Enteric glia 1', 'Enteric glia 2', 'Enteric glia 3', 'Enteric glia 4', 'Enteric glia 5', 'Enteric glia 6', 'Enteric glia 7', 'Enteric glia 8', 'Enteric neurons 1', 'Enteric neurons 2', 'Enteric neurons 3', 'Endothelial 1', 'Endothelial 2', 'Vascular endothelial 1', 'Vascular endothelial 2', 'Vascular endothelial 3', 'Vascular endothelial 4', 'Cajal-Retzius cells', 'Purkinje cells', 'Interneurons 1', 'Interneurons 2', 'Interneurons 3', 'Glutamatergic neurons 1', 'Glutamatergic neurons 2', 'Glutamatergic neurons 3', 'Glutamatergic neurons 4', 'Glutamatergic neurons 5', 'Glutamatergic neurons 6', 'Neurons 1', 'Neurons 2', 'Hypoth. inhib. neurons', 'Hypoth. neurons', 'Intermediate progenitors 1', 'Intermediate progenitors 2', 'Neuroblasts 1', 'Neuroblasts 2', 'Neuroblasts 3', 'Neuroblasts 4', 'Neuroblasts 5', 'Neuroblasts 6', 'Neuroblasts 7', 'Neuroblasts 8')
cell_type_order_karolinska_lvl1 <- c('Pyramidal (CA1)','Pyramidal (SS)','Interneurons','Medium spiny neuron', 'Astrocytes / Ependymal','Dopaminergic neuron','DA neuroblast','Embr. DA neuron','Embr. GABA neuron','Embr. midbrain nucl. neuron','Endothelial-Mural','Microglia','Neural progenitors','Neuroblasts','Oligodendrocyte precursors','Oligodendrocytes','Radial glia like','Striatal interneuron','Vasc. leptomeningeal cells', 'Hypoth. DA neurons','Hypoth. GABAergic neuron','Hypoth. glutamat. neuron', 'Oxytocin/vasopressin neurons','Serotonergic neuron')
cell_type_order_macparland <- c('Antibody secreting B cells','CD3+ αβ T cells','Central venous LSECs','Cholangiocytes','Erthyroid cells','Hep 1','Hep 2','Hep 3','Hep 4','Hep 5','Hep 6','Inflammatory macrophages','Mature B cells','NK-like cells','Non-inflammatory Macrophages','Periportal LSECs','Portal endothelial cells','Stellate cells','γδ T cells 1','γδ T cells 2')

file_path_10x_genomics_lvl1 <- '~/umcu_internship/Single-Cell/data/10x_Genomics_16_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_10x_genomics_lvl2 <- '~/umcu_internship/Single-Cell/data/10x_Genomics_59_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_karolinska_lvl1 <- '~/umcu_internship/Single-Cell/data/Karolinska_level_1_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_karolinska_lvl2 <- '~/umcu_internship/Single-Cell/data/Karolinska_level_2_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_macparland <- '~/umcu_internship/Single-Cell/data/MacParland_20_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'


########## Plot DEPICT results for 10x Genomics, k-means clustering (k=20), 16 cell types ##########
# (1) Load data
heatmap_10x_genomics <- read.csv(file_path_10x_genomics_lvl1)
# (2) Mark significant P-values (Bonferroni-threshold: 0.05/(16*6) )
significant_frames_10x = heatmap_10x_genomics[heatmap_10x_genomics$Significant, c('GWAS', 'Cell.type')]
significant_frames_10x <- significant_frames_10x %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(Cell.type = as.character(Cell.type)) %>%     
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Neuroblasts (2)', which(cell_type_order_10x_genomics_lvl1 == 'Neuroblasts (2)')))
significant_frames_10x <- data.frame(lapply(significant_frames_10x,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_10x_genomics) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_10x_genomics_lvl1), fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_10x, size = 0.75, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(14) +
    theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
    theme(title = element_text(size = 9), axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for 10x Genomics, graph-based clustering, 59 cell types ##########
# (1) Load data
heatmap_10x_genomics <- read.csv(file_path_10x_genomics_lvl2)
# (2) Mark significant P-values (Bonferroni-threshold: 0.05/(59*6) )
significant_frames_10x = heatmap_10x_genomics[heatmap_10x_genomics$Significant, c('GWAS', 'Cell.type')]
significant_frames_10x <- significant_frames_10x %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(Cell.type = as.character(Cell.type)) %>%     
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Hypoth. neurons', which(cell_type_order_10x_genomics_lvl2 == 'Hypoth. neurons'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Hypoth. inhib. neurons', which(cell_type_order_10x_genomics_lvl2 == 'Hypoth. inhib. neurons'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Neuroblasts 5', which(cell_type_order_10x_genomics_lvl2 == 'Neuroblasts 5'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Radial glia 1', which(cell_type_order_10x_genomics_lvl2 == 'Radial glia 1'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Glutamatergic neurons 5', which(cell_type_order_10x_genomics_lvl2 == 'Glutamatergic neurons 5')))
significant_frames_10x <- data.frame(lapply(significant_frames_10x,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_10x_genomics) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_10x_genomics_lvl2), fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_10x, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(title = element_text(size = 9), axis.text.y = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for Karolinska Institute dataset, 24 cell types ##########
# (1) Load data
heatmap_karolinska <- read.csv(file_path_karolinska_lvl1)
# (2) Mark significant P-values (Bonferroni-threshold: 0.05/(24*6) )
significant_frames_karolinska = heatmap_karolinska[heatmap_karolinska$Significant, c('GWAS', 'Cell.type')]
significant_frames_karolinska <- significant_frames_karolinska %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'Breast cancer (2017)', which(gwas_order == 'Breast cancer (2017)'))) %>% 
    mutate(Cell.type = as.character(Cell.type)) %>%     
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Interneurons', which(cell_type_order_karolinska_lvl1 == 'Interneurons'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Medium spiny neuron', which(cell_type_order_karolinska_lvl1 == 'Medium spiny neuron'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Pyramidal (CA1)', which(cell_type_order_karolinska_lvl1 == 'Pyramidal (CA1)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Pyramidal (SS)', which(cell_type_order_karolinska_lvl1 == 'Pyramidal (SS)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Radial glia like', which(cell_type_order_karolinska_lvl1 == 'Radial glia like'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Neural progenitors', which(cell_type_order_karolinska_lvl1 == 'Neural progenitors')))
significant_frames_karolinska <- data.frame(lapply(significant_frames_karolinska,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_karolinska) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_karolinska_lvl1), fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_karolinska, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for Karolinska Institute dataset, 149 cell types ##########
# (1) Load data
heatmap_karolinska <- read.csv(file_path_karolinska_lvl2)
# (2) Mark significant P-values (Bonferroni-threshold: 0.05/(149*6) )
significant_frames_karolinska = heatmap_karolinska[heatmap_karolinska$Significant, c('GWAS', 'Cell.type')]
significant_frames_karolinska <- significant_frames_karolinska %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)')))
significant_frames_karolinska <- data.frame(lapply(significant_frames_karolinska,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_karolinska) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = Cell.type, fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_karolinska, size = 0.25, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(title = 'KI dataset (level 2)', x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.y = element_text(size=5), axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for MacParland human liver dataset, 20 cell types ##########
# (1) Load data
heatmap_macparland <- read.csv(file_path_macparland)
# (2) Mark significant P-values (Bonferroni-threshold: 0.05/(20*6) )
significant_frames_macparland = heatmap_macparland[heatmap_macparland$Significant, c('GWAS', 'Cell.type')]
significant_frames_macparland$GWAS = which(gwas_order == 'Intelligence (2018)')
significant_frames_macparland$Cell.type = which(cell_type_order_macparland == 'γδ T cells 2')
# (3) Plot heatmap
ggplot(data=heatmap_macparland) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_macparland), fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_macparland, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))
