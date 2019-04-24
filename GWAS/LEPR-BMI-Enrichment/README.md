![image](https://user-images.githubusercontent.com/24732704/55021982-f822ec00-4ff9-11e9-802a-649cfdb4892c.png)

## UMC Utrecht internship code repository - GWAS - Enrichment of specific genes in summary statistics

### Background
Script used to determine whether a given set of genes that are enriched in cells expressing LEPR are also enriched in BMI GWAS summary statistics.

### Scripts
* **LEPR_genes_BMI_GWAS_enrichment.R** - Script to map mouse genes to 1:1 human orthologs and compare human genes against MAGMA gene analysis output for BMI GWAS.

### Requirements
- Unix environment (I used Ubuntu 18.04)
- R version 3.5> (https://www.r-project.org/) with the [One2One](https://github.com/NathanSkene/One2One), [readxl](https://cran.r-project.org/web/packages/readxl/) and [xlsx](https://cran.r-project.org/web/packages/xlsx/) packages installed
- File with list of mouse genes of interest (specifically an Excel sheet, yet this can be changed depending on data format)
- MAGMA gene analysis output (magma.genes.out)

### Usage
**mouse_to_human_gene_mapping.R**:
1. `Rscript LEPR_genes_BMI_GWAS_enrichment.R`
