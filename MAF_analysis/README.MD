![image](https://user-images.githubusercontent.com/24732704/55021982-f822ec00-4ff9-11e9-802a-649cfdb4892c.png)

## UMC Utrecht internship code repository - Analysis of MDD MAF

### Background
MAF indicates how common a SNP is in a population, and can be used to discern rare- from common variants. While processing summary statistics for the 2018 MDD GWAS [1], I noticed MAFs were given for both case- and control populations. Before deciding on how to deal with this, and what filter to apply, I wanted to check the distribution of MAFs in MDD and control populations.

### Scripts
- **maf_comparison.sh** - Script that downloads the summary statistics, unzips it, and selects specific columns of data
- **visualize_mdd_maf.R** - Script that loads the data and visualizes the distribution of MAFs for case and control populations in a boxplot and violin plot

### Data
- MDD 2018 GWAS summary statistics, will automatically be downloaded and processed by **maf_comparison.sh**

### Usage

1. `chmod 755 maf_comparison.sh`
1. `./maf_comparison.sh`
1. Open **visualize_mdd_maf.R** with RStudio
1. Set FILE_PATH to the location of the output on your system
1. Run the code (loading data might take a moment!) and export the plots as you please.

### Requirements
- Unix environment (I used Ubuntu 18.04) with wget, gunzip & awk installed
- R version 3.5> (https://www.r-project.org/) with the ggplot2 package
- RStudio version 1.1> (https://www.rstudio.com/)

### References
[1] Wray et al. (2018), "Genome-wide association analyses identify 44 risk variants and refine the genetic architecture of major depression", _Nature Genetics_ 50, p668-681, doi: https://doi.org/10.1038/s41588-018-0090-3

### Glossary
* **MAF** - Minor Allele Frequency; the frequency at which the second most common allele occur in a population
* **MDD** - Major Depressive Disorder; mental disorder characterized by low mood, low self-esteem and loss of interest in normal activities
* **SNP** - Single base-pair change that occurs on a specific genome position
