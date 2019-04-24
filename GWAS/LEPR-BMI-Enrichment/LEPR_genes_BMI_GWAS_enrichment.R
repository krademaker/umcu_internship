# TITLE:    LEPR_genes_BMI_GWAS_enrichment.R
# ABOUT:    Script to (1) map mouse genes to 1:1 human orthologs & (2) compare human genes against MAGMA gene-
#           analysis output for BMI GWAS.
# INPUT:    IPIN_0.01_up.xlsx - Excel file containing a single column with mouse genes in the sheet 'Gene list'.
# INPUT:    magma.genes.out - MAGMA gene analysis output for Yengo et al. (2018) BMI GWAS.
# OUTPUT:   LEPR_genes_BMI_GWAS_enrichment.xlsx - File containing sheets with (1) MAGMA gene analysis output-
#           for selected genes & (2) comparison of mouse and human gene names for 1:1 orthologs.
# AUTHOR:   Koen Rademaker
# DATE:     24 April 2019

########## Load required libraries ##########
library(One2One)
library(readxl)
library(xlsx)

########## Read list of mouse genes enriched in LEPR-expressing cells ##########
mouse_gene_list <- read_excel('IPIN_0.01_up.xlsx',
                              sheet = 'Gene list',
                              col_names = 'Gene name')

########## Read MAGMA gene analysis output for BMI GWAS ##########
bmi_gwas_magma_output <- read.delim('magma.genes.out')

########## Load orthology data for mouse and human ##########
all_homologs = load.homologs()
species_1 = 'human'
species_2 = 'mouse'
orthology_data = analyse.orthology(species_1, species_2, all_homologs)

########## Get 1:1 orthologs between mouse and human, extract list of human genes ##########
intersecting_genes <- intersect(mouse_gene_list$`Gene name`, orthology_data$orthologs_one2one$mouse.symbol)
subset_orthology_data <- filter(orthology_data$orthologs_one2one, mouse.symbol %in% intersecting_genes)

########## Re-arrange and alphabetically sort on mouse gene names ##########
subset_orthology_data <- subset_orthology_data[3:2]
subset_orthology_data <- subset_orthology_data[order(subset_orthology_data$mouse.symbol),]

########## Extract MAGMA output for human ortholog genes ########## 
subset_bmi_gwas_magma_output <- merge(x = bmi_gwas_magma_output, y = subset_orthology_data,
                                      by.x = 'SYMBOL', by.y = 'human.symbol')

########## Write output to Excel file ##########
filename_output <- 'LEPR_genes_BMI_GWAS_enrichment.xlsx'
write.xlsx(x = subset_bmi_gwas_magma_output, file = filename_output,
           sheetName = 'BMI_GWAS_MAGMA_Gene_Analysis', col.names = TRUE,
           row.names = FALSE)
write.xlsx(x = subset_orthology_data, file = filename_output,
           sheetName='Mouse_Human_Genes', col.names = TRUE,
           row.names = FALSE, append = TRUE)
