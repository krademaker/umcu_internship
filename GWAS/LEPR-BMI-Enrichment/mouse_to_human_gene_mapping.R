# TITLE:    mouse_to_human_gene_mapping.R
# ABOUT:    Script to load a list of mouse genes from an input file, and subsequently convert them a list of 1:1 ortholog human genes.
# INPUT:    [filename].xlsx - List of mouse genes, stored as a single column in an Excel file.
# OUTPUT:   [filename].txt - List of human genes.
# AUTHOR:   Koen Rademaker
# DATE:     23 April 2019

########## Load required libraries ##########
library(One2One)
library(readxl)

########## Read input data ##########
mouse_gene_list <- read_excel('~/IPIN_0.01_up.xlsx', sheet = 'Gene list', col_names = 'Gene name')

########## Load orthology for mice and humans ##########
all_homologs = load.homologs()
species_1 = 'human'
species_2 = 'mouse'
orthology_data = analyse.orthology(species_1, species_2, all_homologs)

########## Get 1:1 orthologs and extract human gene names ##########
intersecting_genes <- intersect(mouse_gene_list$`Gene name`, orthology_data$orthologs_one2one$mouse.symbol)
intersecting_genes_annotation <- filter(orthology_data$orthologs_one2one, mouse.symbol %in% intersecting_genes)
human_gene_list <- intersecting_genes_annotation$human.symbol

########## Write human gene list to a text file ##########
write(human_gene_list, file='human_gene_orthologs.txt', sep='\n')
