# TITLE:    plot_combined_gwas_depict_results.R
# ABOUT:    Script to plot DEPICT results and highlight significant cell types for different scRNA-seq datasets.
# INPUT:    file_path_10x_genomics: Combined GWAS DEPICT results for 10x Genomics dataset
# INPUT:    file_path_karolinska: Combined GWAS DEPICT results for Karolinska Institute dataset
# INPUT:    file_path_macparland: Combined GWAS DEPICT results for MacParland dataset
# AUTHOR:   Koen Rademaker
# DATE:     6 June 2019


########## Load required libraries ##########
library(ggplot2)
library(dplyr)


########## Set variables ##########
gwas_order <- c('SCZ (2018)', 'EA (2018)', 'MDD (2018)', 'Intelligence (2018)', 'BMI (2015)', 'Breast cancer (2017)')
karolinska_cell_type_order <- c('Pyramidal (CA1)','Pyramidal (SS)','Interneurons','Medium spiny neuron', 'Astrocytes / Ependymal','Dopaminergic neuron','DA neuroblast','Embr. DA neuron','Embr. GABA neuron','Embr. midbrain nucl. neuron','Endothelial-Mural','Microglia','Neural progenitors','Neuroblasts','Oligodendrocyte precursors','Oligodendrocytes','Radial glia like','Striatal interneuron','Vasc. leptomeningeal cells', 'Hypoth. DA neurons','Hypoth. GABAergic neuron','Hypoth. glutamat. neuron', 'Oxytocin/vasopressin neurons','Serotonergic neuron')
macparland_cell_type_order <- c('Antibody secreting B cells','CD3+ αβ T cells','Central venous LSECs','Cholangiocytes','Erthyroid cells','Hep 1','Hep 2','Hep 3','Hep 4','Hep 5','Hep 6','Inflammatory macrophages','Mature B cells','NK-like cells','Non-inflammatory Macrophages','Periportal LSECs','Portal endothelial cells','Stellate cells','γδ T cells 1','γδ T cells 2')

file_path_10x_genomics <- '~/umcu_internship/Single-Cell/data/10x_Genomics_16_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_karolinska <- '~/umcu_internship/Single-Cell/data/Karolinska_level_1_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_macparland <- '~/umcu_internship/Single-Cell/data/MacParland_20_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'


########## Plot DEPICT results for 10x Genomics, k-means clustering (k=20), 16 cell types ##########
heatmap_10x_genomics <- read.csv(file_path_10x_genomics)

significant_frames_10x = heatmap_10x_genomics[heatmap_10x_genomics$Significant, c('GWAS', 'Cell.type')]
significant_frames_10x$GWAS = which(gwas_order == 'EA (2018)')
significant_frames_10x$Cell.type = as.integer(significant_frames_10x$Cell.type)

ggplot(data=heatmap_10x_genomics) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = Cell.type, fill = P)) +
    scale_fill_gradient2(low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = 'colourbar') +
    geom_rect(data = significant_frames_10x, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(title = 'DEPICT results for 10x Genomics dataset (108K subset)', x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for Karolinska Institute dataset, 24 cell types ##########
heatmap_karolinska <- read.csv(file_path_karolinska)

significant_frames_karolinska = heatmap_karolinska[heatmap_karolinska$Significant, c('GWAS', 'Cell.type')]
significant_frames_karolinska <- significant_frames_karolinska %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'Breast cancer (2017)', which(gwas_order == 'Breast cancer (2017)'))) %>% 
    mutate(Cell.type = as.character(Cell.type)) %>%     
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Interneurons', which(karolinska_cell_type_order == 'Interneurons'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Medium spiny neuron', which(karolinska_cell_type_order == 'Medium spiny neuron'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Pyramidal (CA1)', which(karolinska_cell_type_order == 'Pyramidal (CA1)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Pyramidal (SS)', which(karolinska_cell_type_order == 'Pyramidal (SS)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Radial glia like', which(karolinska_cell_type_order == 'Radial glia like'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Neural progenitors', which(karolinska_cell_type_order == 'Neural progenitors')))
significant_frames_karolinska <- data.frame(lapply(significant_frames_karolinska,as.numeric))

ggplot(data=heatmap_karolinska) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = karolinska_cell_type_order), fill = P)) +
    scale_fill_gradient2(low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = 'colourbar') +
    geom_rect(data = significant_frames_karolinska, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(title = 'DEPICT results for Karolinska Institute dataset', x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for MacParland human liver dataset, 20 cell types ##########
heatmap_macparland <- read.csv(file_path_macparland)

significant_frames_macparland = heatmap_macparland[heatmap_macparland$Significant, c('GWAS', 'Cell.type')]
significant_frames_macparland$GWAS = which(gwas_order == 'Intelligence (2018)')
significant_frames_macparland$Cell.type = which(macparland_cell_type_order == 'γδ T cells 2')

ggplot(data=heatmap_macparland) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = macparland_cell_type_order), fill = P)) +
    scale_fill_gradient2(low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = 'colourbar') +
    geom_rect(data = significant_frames_macparland, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(title = 'DEPICT results for human liver dataset', x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))