# TITLE:    plot_combined_gwas_magma_results.R
# ABOUT:    Script to plot MAGMA results and highlight significant cell types for different scRNA-seq datasets.
# INPUT:    file_path_10x_genomics_linear: Combined GWAS MAGMA linear mode results for 10x Genomics dataset (16 cell types)
# INPUT:    file_path_10x_genomics_top10: Combined GWAS MAGMA top 10% mode results for 10x Genomics dataset (16 cell types)
# INPUT:    file_path_karolinska_linear: Combined GWAS MAGMA linear mode results for KI dataset level 1 (24 cell types)
# INPUT:    file_path_karolinska_top10: Combined GWAS MAGMA top 10% mode results for KI dataset level 1 (24 cell types)
# INPUT:    file_path_macparland_linear: Combined GWAS MAGMA linear mode results for MacParland dataset (20 cell types)
# INPUT:    file_path_macparland_top10: Combined GWAS MAGMA top 10% mode results for MacParland dataset (20 cell types)
# AUTHOR:   Koen Rademaker
# DATE:     20 June 2019


########## Load required libraries ##########
library(ggplot2)
library(dplyr)


########## Set variables ##########
gwas_order <- c('SCZ (2018)', 'EA (2018)', 'MDD (2018)', 'Intelligence (2018)', 'BMI (2015)', 'Breast cancer (2017)')
cell_type_order_10x_genomics <- c('Astrocytes 1', 'Astrocytes 2', 'Cajal-Retzius cells', 'Endothelial cells', 'Enteric glial cells', 'Enteric neurons', 'Glutamatergic neurons', 'Intermediate progenitors', 'Interneurons 1', 'Interneurons 2', 'Microglia', 'Neuroblasts 1', 'Neuroblasts 2', 'Neurons', 'Oligodendrocytes', 'Vascular endothelial cells')
cell_type_order_karolinska_lvl1 <- c('Pyramidal (CA1)','Pyramidal (SS)','Interneurons','Medium spiny neuron', 'Astrocytes / Ependymal','Dopaminergic neuron','DA neuroblast','Embr. DA neuron','Embr. GABA neuron','Embr. midbrain nucl. neuron','Endothelial-Mural','Microglia','Neural progenitors','Neuroblasts','Oligodendrocyte precursors','Oligodendrocytes','Radial glia like','Striatal interneuron','Vasc. leptomeningeal cells', 'Hypoth. DA neurons','Hypoth. GABAergic neuron','Hypoth. glutamat. neuron', 'Oxytocin/vasopressin neurons','Serotonergic neuron')
cell_type_order_macparland <- c('Antibody secreting B cells','CD3+ αβ T cells','Central venous LSECs','Cholangiocytes','Erthyroid cells','Hep 1','Hep 2','Hep 3','Hep 4','Hep 5','Hep 6','Inflammatory macrophages','Mature B cells','NK-like cells','Non-inflammatory Macrophages','Periportal LSECs','Portal endothelial cells','Stellate cells','γδ T cells 1','γδ T cells 2')

file_path_10x_genomics_linear <- '~/umcu_internship/Single-Cell/data/10x_Genomics_combined_GWAS_MAGMA_linear_results_heatmap.csv'
file_path_10x_genomics_top10 <- '~/umcu_internship/Single-Cell/data/10x_Genomics_combined_GWAS_MAGMA_top10_results_heatmap.csv'
file_path_karolinska_linear <- '~/umcu_internship/Single-Cell/data/KI_level_1_combined_GWAS_MAGMA_linear_results_heatmap.csv'
file_path_karolinska_top10 <- '~/umcu_internship/Single-Cell/data/KI_level_1_combined_GWAS_MAGMA_top10_results_heatmap.csv'
file_path_macparland_linear <- '~/umcu_internship/Single-Cell/data/MacParland_combined_GWAS_MAGMA_linear_results_heatmap.csv'
file_path_macparland_top10 <- '~/umcu_internship/Single-Cell/data/MacParland_combined_GWAS_MAGMA_top10_results_heatmap.csv'


########## Plot MAGMA linear mode results for 10x Genomics dataset ##########
# (1) Load data
heatmap_10x_linear <- read.csv(file_path_10x_genomics_linear)
# (2) Mark significant P-values (Bonferroni-threshold: 0.05/(16*6) )
significant_frames_10x_linear = heatmap_10x_linear[heatmap_10x_linear$Significant, c('GWAS', 'Cell.type')]
significant_frames_10x_linear <- significant_frames_10x_linear %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'MDD (2018)', which(gwas_order == 'MDD (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'Intelligence (2018)', which(gwas_order == 'Intelligence (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'BMI (2015)', which(gwas_order == 'BMI (2015)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'Breast cancer (2017)', which(gwas_order == 'Breast cancer (2017)')))
significant_frames_10x_linear$Cell.type = as.integer(significant_frames_10x_linear$Cell.type)
significant_frames_10x_linear <- data.frame(lapply(significant_frames_10x_linear,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_10x_linear) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = Cell.type, fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_10x_linear, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(title = element_text(size = 9), axis.text.y = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot MAGMA linear mode results for KI dataset ##########
# (1) Load data
heatmap_ki_linear <- read.csv(file_path_karolinska_linear)
# (2) Mark significant P-values (Bonferroni-threshold: 0.05/(24*6) )
significant_frames_ki_linear = heatmap_ki_linear[heatmap_ki_linear$Significant, c('GWAS', 'Cell.type')]
significant_frames_ki_linear <- significant_frames_ki_linear %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'BMI (2015)', which(gwas_order == 'BMI (2015)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'Intelligence (2018)', which(gwas_order == 'Intelligence (2018)')))  %>%
    mutate(Cell.type = as.character(Cell.type)) %>%     
    mutate(Cell.type = replace(Cell.type, Cell.type == 'DA neuroblast', which(cell_type_order_karolinska_lvl1 == 'DA neuroblast'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Embr. GABA neuron', which(cell_type_order_karolinska_lvl1 == 'Embr. GABA neuron'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Embr. midbrain nucl. neuron', which(cell_type_order_karolinska_lvl1 == 'Embr. midbrain nucl. neuron'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Interneurons', which(cell_type_order_karolinska_lvl1 == 'Interneurons'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Medium spiny neuron', which(cell_type_order_karolinska_lvl1 == 'Medium spiny neuron'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Pyramidal (CA1)', which(cell_type_order_karolinska_lvl1 == 'Pyramidal (CA1)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Pyramidal (SS)', which(cell_type_order_karolinska_lvl1 == 'Pyramidal (SS)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Serotonergic neuron', which(cell_type_order_karolinska_lvl1 == 'Serotonergic neuron'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Striatal interneuron', which(cell_type_order_karolinska_lvl1 == 'Striatal interneuron'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Embr. DA neuron', which(cell_type_order_karolinska_lvl1 == 'Embr. DA neuron'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Neuroblasts', which(cell_type_order_karolinska_lvl1 == 'Neuroblasts'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Hypoth. DA neurons', which(cell_type_order_karolinska_lvl1 == 'Hypoth. DA neurons'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Oligodendrocytes', which(cell_type_order_karolinska_lvl1 == 'Oligodendrocytes')))
significant_frames_ki_linear <- data.frame(lapply(significant_frames_ki_linear,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_ki_linear) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_karolinska_lvl1), fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_ki_linear, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot MAGMA linear mode results for MacParland dataset ##########
# (1) Load data
heatmap_macparland_linear <- read.csv(file_path_macparland_linear)
# (2) Mark significant P-values (Bonferroni-threshold: 0.05 / (20*6) )
significant_frames_macparland_linear = heatmap_macparland_linear[heatmap_macparland_linear$Significant, c('GWAS', 'Cell.type')]
significant_frames_macparland_linear <- significant_frames_macparland_linear %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'Intelligence (2018)', which(gwas_order == 'Intelligence (2018)')))  %>%
    mutate(Cell.type = as.character(Cell.type)) %>%     
    mutate(Cell.type = replace(Cell.type, Cell.type == 'γδ T cells 2', which(cell_type_order_macparland == 'γδ T cells 2'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Erthyroid cells', which(cell_type_order_macparland == 'Erthyroid cells'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Central venous LSECs', which(cell_type_order_macparland == 'Central venous LSECs'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'CD3+ αβ T cells', which(cell_type_order_macparland == 'CD3+ αβ T cells'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'NK-like cells', which(cell_type_order_macparland == 'NK-like cells'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'γδ T cells 1', which(cell_type_order_macparland == 'γδ T cells 1')))
significant_frames_macparland_linear <- data.frame(lapply(significant_frames_macparland_linear,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_macparland_linear) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_macparland), fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'brown4', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_macparland_linear, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))

########## Plot MAGMA top 10% mode results for 10x Genomics dataset ##########
# (1) Load data
heatmap_10x_top10 <- read.csv(file_path_10x_genomics_top10)
# (2) Mark significant P-values (Bonferroni threshold: 0.05/(16*6) )
significant_frames_10x_top10 = heatmap_10x_top10[heatmap_10x_top10$Significant, c('GWAS', 'Cell.type')]
significant_frames_10x_top10 <- significant_frames_10x_top10 %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'Intelligence (2018)', which(gwas_order == 'Intelligence (2018)')))
significant_frames_10x_top10$Cell.type = as.integer(significant_frames_10x_top10$Cell.type)
significant_frames_10x_top10 <- data.frame(lapply(significant_frames_10x_top10,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_10x_top10) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = Cell.type, fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_10x_top10, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(title = element_text(size = 9), axis.text.y = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot MAGMA top 10% mode results for KI dataset ##########
# (1) Load data
heatmap_ki_top10 <- read.csv(file_path_karolinska_top10)
# (2) Mark significant P-values (Bonferroni threshold: 0.05/(24*6) )
significant_frames_ki_top10 <- heatmap_ki_top10[heatmap_ki_top10$Significant, c('GWAS', 'Cell.type')]
significant_frames_ki_top10 <- significant_frames_ki_top10 %>% 
    mutate(GWAS = as.character(GWAS)) %>% 
    mutate(GWAS = replace(GWAS, GWAS == 'SCZ (2018)', which(gwas_order == 'SCZ (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'EA (2018)', which(gwas_order == 'EA (2018)'))) %>%
    mutate(GWAS = replace(GWAS, GWAS == 'Intelligence (2018)', which(gwas_order == 'Intelligence (2018)'))) %>%
    mutate(Cell.type = as.character(Cell.type)) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Medium spiny neuron', which(cell_type_order_karolinska_lvl1 == 'Medium spiny neuron'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Pyramidal (CA1)', which(cell_type_order_karolinska_lvl1 == 'Pyramidal (CA1)'))) %>%
    mutate(Cell.type = replace(Cell.type, Cell.type == 'Pyramidal (SS)', which(cell_type_order_karolinska_lvl1 == 'Pyramidal (SS)')))
significant_frames_ki_top10 <- data.frame(lapply(significant_frames_ki_top10,as.numeric))
# (3) Plot heatmap
ggplot(data=heatmap_ki_top10) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_karolinska_lvl1), fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'navyblue', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    geom_rect(data = significant_frames_ki_top10, size = 0.5, fill = NA, colour = 'red', aes(xmin = GWAS - 0.5, xmax = GWAS + 0.5, ymin = Cell.type - 0.5, ymax = Cell.type + 0.5)) +
    labs(x = 'GWAS', y = 'Cell type') +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


########## Plot MAGMA top 10% mode results for MacParland dataset ##########
# (1) Load data
heatmap_macparland_top10 <- read.csv(file_path_macparland_top10)
# (2) Plot heatmap
ggplot(data=heatmap_macparland_top10) +
    geom_raster(aes(x = factor(GWAS, level = gwas_order), y = factor(Cell.type, level = cell_type_order_macparland), fill = P)) +
    scale_fill_gradient2(limits = c(0.0, 1.0), low = 'brown4', mid = 'snow', high = 'white', midpoint = 0.1, space = 'Lab', na.value = 'grey50', guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
    labs(x = 'GWAS', y = 'Cell type')
    theme_classic(12) +
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))