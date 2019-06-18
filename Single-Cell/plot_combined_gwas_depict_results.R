# TITLE:    plot_combined_gwas_depict_results.R
# ABOUT:    Script to plot DEPICT results and highlight significant cell types for different scRNA-seq datasets.
# INPUT:    file_path_10x_genomics_lvl1: Combined GWAS DEPICT results for 10x Genomics dataset level 1 (n=16)
# INPUT:    file_path_10x_genomics_lvl2: Combined GWAS DEPICT results for 10x Genomicd dataset level 2 (n=59)
# INPUT:    file_path_karolinska_lvl1: Combined GWAS DEPICT results for Karolinska Institute (KI) dataset level 1 (n=24)
# INPUT:    file_path_karolinska_lvl2: Combined GWAS DEPICT results for Karolinska Institute (KI) dataset levle 2 (n=149)
# INPUT:    file_path_macparland: Combined GWAS DEPICT results for MacParland dataset
# AUTHOR:   Koen Rademaker
# DATE:     18 June 2019


########## Load required libraries ##########
library(ggplot2)
library(dplyr)


########## Set variables ##########
gwas_order <- c('SCZ (2018)', 'EA (2018)', 'MDD (2018)', 'Intelligence (2018)', 'BMI (2015)', 'Breast cancer (2017)')
cell_type_order_10x_genomics <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59')
cell_type_order_karolinska_lvl1 <- c('Pyramidal (CA1)','Pyramidal (SS)','Interneurons','Medium spiny neuron', 'Astrocytes / Ependymal','Dopaminergic neuron','DA neuroblast','Embr. DA neuron','Embr. GABA neuron','Embr. midbrain nucl. neuron','Endothelial-Mural','Microglia','Neural progenitors','Neuroblasts','Oligodendrocyte precursors','Oligodendrocytes','Radial glia like','Striatal interneuron','Vasc. leptomeningeal cells', 'Hypoth. DA neurons','Hypoth. GABAergic neuron','Hypoth. glutamat. neuron', 'Oxytocin/vasopressin neurons','Serotonergic neuron')
cell_type_order_macparland <- c('Antibody secreting B cells','CD3+ αβ T cells','Central venous LSECs','Cholangiocytes','Erthyroid cells','Hep 1','Hep 2','Hep 3','Hep 4','Hep 5','Hep 6','Inflammatory macrophages','Mature B cells','NK-like cells','Non-inflammatory Macrophages','Periportal LSECs','Portal endothelial cells','Stellate cells','γδ T cells 1','γδ T cells 2')

file_path_10x_genomics_lvl1 <- '~/umcu_internship/Single-Cell/data/10x_Genomics_16_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_10x_genomics_lvl2 <- '~/umcu_internship/Single-Cell/data/10x_Genomics_59_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_karolinska_lvl1 <- '~/umcu_internship/Single-Cell/data/Karolinska_level_1_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_karolinska_lvl2 <- '~/umcu_internship/Single-Cell/data/Karolinska_level_2_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'
file_path_macparland <- '~/umcu_internship/Single-Cell/data/MacParland_20_cell_types_combined_GWAS_DEPICT_results_heatmap.csv'


########## Plot DEPICT results for 10x Genomics, k-means clustering (k=20), 16 cell types ##########
heatmap_10x_genomics <- read.csv(file_path_10x_genomics_lvl1)

significant_frames_10x = heatmap_10x_genomics[heatmap_10x_genomics$Significant, c('GWAS', 'Cell.type')]
significant_frames_10x$GWAS = which(gwas_order == 'EA (2018)')
significant_frames_10x$Cell.type = as.integer(significant_frames_10x$Cell.type)

ggplot(data=heatmap_10x_genomics) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = Cell.type, fill = P)) +
    scale_fill_gradient2(low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = 'colourbar') +
    geom_rect(data = significant_frames_10x, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(title = 'DEPICT results for 10x Genomics dataset (k-means clustering)', x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(title = element_text(size = 9), axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for 10x Genomics, graph-based clustering, 59 cell types ##########
heatmap_10x_genomics <- read.csv(file_path_10x_genomics_lvl2)

significant_frames_10x = heatmap_10x_genomics[heatmap_10x_genomics$Significant, c('GWAS', 'Cell.type')]
significant_frames_10x <- significant_frames_10x %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(Cell.type = as.character(Cell.type)) %>%     
    mutate(Cell.type = replace(Cell.type, Cell.type == '19', which(cell_type_order_10x_genomics == '19'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == '20', which(cell_type_order_10x_genomics == '20'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == '25', which(cell_type_order_10x_genomics == '25'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == '30', which(cell_type_order_10x_genomics == '30'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == '43', which(cell_type_order_10x_genomics == '43')))
significant_frames_10x <- data.frame(lapply(significant_frames_10x,as.numeric))

ggplot(data=heatmap_10x_genomics) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_10x_genomics), fill = P)) +
    scale_fill_gradient2(low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = 'colourbar') +
    geom_rect(data = significant_frames_10x, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(title = 'DEPICT results for 10x Genomics dataset (graph-based clustering)', x = 'GWAS', y = 'Cluster') +
    theme_classic(12) +
    theme(title = element_text(size = 9), axis.text.y = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for Karolinska Institute dataset, 24 cell types ##########
heatmap_karolinska <- read.csv(file_path_karolinska_lvl1)

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

ggplot(data=heatmap_karolinska) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_karolinska_lvl1), fill = P)) +
    scale_fill_gradient2(low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = 'colourbar') +
    geom_rect(data = significant_frames_karolinska, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(title = 'KI dataset (level 1)', x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for Karolinska Institute dataset, 149 cell types ##########
heatmap_karolinska <- read.csv(file_path_karolinska_lvl2)

significant_frames_karolinska = heatmap_karolinska[heatmap_karolinska$Significant, c('GWAS', 'Cell.type')]
significant_frames_karolinska <- significant_frames_karolinska %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)')))
significant_frames_karolinska <- data.frame(lapply(significant_frames_karolinska,as.numeric))

ggplot(data=heatmap_karolinska) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = Cell.type, fill = P)) +
    scale_fill_gradient2(low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = 'colourbar') +
    geom_rect(data = significant_frames_karolinska, size = 0.25, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(title = 'KI dataset (level 2)', x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.y = element_text(size=3), axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot DEPICT results for MacParland human liver dataset, 20 cell types ##########
heatmap_macparland <- read.csv(file_path_macparland)

significant_frames_macparland = heatmap_macparland[heatmap_macparland$Significant, c('GWAS', 'Cell.type')]
significant_frames_macparland$GWAS = which(gwas_order == 'Intelligence (2018)')
significant_frames_macparland$Cell.type = which(cell_type_order_macparland == 'γδ T cells 2')

ggplot(data=heatmap_macparland) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_macparland), fill = P)) +
    scale_fill_gradient2(low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = 'colourbar') +
    geom_rect(data = significant_frames_macparland, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(title = 'MacParland dataset', x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))