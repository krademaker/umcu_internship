![image](https://user-images.githubusercontent.com/24732704/55021982-f822ec00-4ff9-11e9-802a-649cfdb4892c.png)

## UMC Utrecht internship code repository - PLINK

### Background
Small exercise to run a logistic regression model GWAS using PLINK [1], on data from the 1000 Genomes Project [2], where sex is the researched trait.

### Approach
1. Build necessary files, and restructure them to serve as input for PLINK.
2. Run logistic regression model, excluding non-autosomal chromosomes and filter MAF > 0.01
3. Plot results in a Manhattan and QQ plot

### Data
Data from the 1000 Genomes Project [2] phase 1 was obtained from their [website](http://www.internationalgenome.org/data/), this included the .bed, .bam and .fam files.

### Scripts
* **run_plink.sh** - Script that perform the necessary preprocessing and runs the logistic regression model in PLINK.
* **plot_plink.sh** - Script that plots the PLINK output in a Manhattan and QQ plot.

### Requirements
- Unix environment (I used Ubuntu 18.04) with awk installed: `sudo apt-get install gawk`
- R version 3.4> (https://www.r-project.org/) with the qqman package
- PLINK version 1.9 (https://www.cog-genomics.org/plink2)

### Usage
- Make sure PLINK is added to your PATH variable so it can be directly called
- `bash run_plink.sh`
- Run the R script, make sure that `sum_stats_path` is set to the correct path.

### References
[1] Purcell et al. (2007), "PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses", _American Journal of Human Genetics_ 81, p559-575, doi: https://doi.org/10.1086/519795 

[2] The 1000 Genomes Project Consortium (2015), "A global reference for human genetic variation", _Nature_ 526, p67-74, doi: https://doi.org/10.1038/nature15393
