![image](https://user-images.githubusercontent.com/24732704/55021982-f822ec00-4ff9-11e9-802a-649cfdb4892c.png)

## UMC Utrecht internship code repository - GWAS - Enrichment of specific genes in summary statistics

### Background
Scripts used to determine whether a given set of genes that are enriched in cells expressing LEPR are also enriched in BMI GWAS summary statistics.

### Scripts
* **mouse_to_human_gene_mapping.R** - Script that converts mouse genes to 1:1 human orthologs and outputs it to a file for later comparison.

### Requirements
- Unix environment (I used Ubuntu 18.04)
- R version 3.5> (https://www.r-project.org/) with the [One2One](https://github.com/NathanSkene/One2One) and [readxl](https://cran.r-project.org/web/packages/readxl/) packages installed
- File with list of mouse genes of interest (specifically an Excel sheet, yet this can be changed depending on data format)

### Usage
**mouse_to_human_gene_mapping.R**:
1. `Rscript mouse_to_human_gene_mapping.R`

### Glossary
* **GWAS** - _Genome Wide Association Study_; hypothesis-free method to identify the genetic variants at whole-genome wide level associated with traits such as diseases
