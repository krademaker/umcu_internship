![image](https://user-images.githubusercontent.com/24732704/55021982-f822ec00-4ff9-11e9-802a-649cfdb4892c.png)

## UMC Utrecht internship code repository - GWAS

### Background
GWASs provide insight into the genetic background of traits and diseases, and results are commonly published in summary statistics. This folder and its' subfolders include scripts related to obtaining and processing summary statistics data.

### Scripts
* **download_summary_statistics.sh** - Script for bulk downloading, unzipping and structuring specific summary statistics files I used for my research.

_(details of subfolders are explained in their own READMEs)_

### Usage
**download_summary_statistics.sh**:
1. `chmod 755 download_summary_statistics.sh`
1. `./download_summary_statistics.sh`

### Requirements
- Unix environment (I used Ubuntu 18.04) with awk installed:
`sudo apt-get install gawk`

### Glossary
* **GWAS** - _Genome Wide Association Study_; hypothesis-free method to identify the genetic variants at whole-genome wide level associated with traits such as diseases
