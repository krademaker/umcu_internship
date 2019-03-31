![image](https://user-images.githubusercontent.com/24732704/55021982-f822ec00-4ff9-11e9-802a-649cfdb4892c.png)

## UMC Utrecht internship code repository - Summary statistics QC

### Background
Summary statistics from GWASs can be insightful for further analyses in e.g. DEPICT, yet quality control (QC) is often needed to reduce the data size and filter out non-relevant associations. Various QC steps, written in `awk`, are combined in this script. If multiple approaches are possible for a specific step, they are all included.

### Scripts
- **sum_stats_qc.sh** - Script with general QC functions for summary statistics that can be adapted to your specific data.

### Data
- Summary statistics (optional: list of high-quality imputed SNPs), download links are typically provided in GWAS publications.
- SNP synonyms file filtered to only include SNPs from the 1000 Genomes Project, files available separately from [MAGMA](https://ctg.cncr.nl/software/magma).

### Usage
1. Adapt the script to your data format and filter requirements. Additionally, intermediate output can be written to files.
1. `chmod 755 sum_stats_qc.sh`
1. `./sum_stats_qc.sh`

### Requirements
- Unix environment (I used Ubuntu 18.04) with awk installed:
`sudo apt-get install gawk`
