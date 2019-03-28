![image](https://user-images.githubusercontent.com/24732704/55021982-f822ec00-4ff9-11e9-802a-649cfdb4892c.png)

## UMC Utrecht internship code repository - Annotate rs IDs for SNPs

### Background
GWAS summary statistics usually include rs IDs for SNPs in the 'SNP' column, however, you might need to annotate them in some cases. This script uses a SNP annotation table obtained from UCSC Table Browser [1] (see _UCSC Data_) to annotate rs IDs to SNPs for summary statistics. An additional issue I encountered was for the 2018 schizophrenia GWAS [2], where the 'SNP' column not only includes the genomic position, but also allele information. To this end, this script temporarily removes allele information, annotates rs IDs, and then re-adds allele information.

### Scripts
- **snp_annotate_rs_ids.sh** - Script that annotates rs IDs for SNPs

### UCSC Data
The following parameters were used to obtain the SNP annotation table from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables):

| Parameter          | Value(s)                                      |
|:-------------------|:----------------------------------------------|
| Clade              | Mammal                                        |
| Genome             | Human                                         |
| Assembly           | Feb. 2009 (GRCh37/hg19)                       |
| Group              | Variation                                     |
| Track              | All SNPs (151)                                |
| Region             | Genome                                        |
| Output format      | Select fields from primary and related tables |
| File type returned | Plain text                                    |
| Selected fields     | chrom, chromStart, chromEnd, name             |

Note: This file easily exceeds 12 GB!

### Data
- UCSC SNP annotation table
- GWAS summary statistics in similar format to that of Pardiñas et al. (2018)

### Usage

1. `chmod 755 snp_annotate_rs_ids.sh`
1. Set file paths in the script
1. `./snp_annotate_rs_ids.sh`
1. Output can be appended to your original summary statistics as you please

### Requirements
- Unix environment (I used Ubuntu 18.04), including AWK:
`sudo apt-get install gawk`

### References
[1] Karolchik et al. (2004), "The UCSC Table Browser data retrieval tool", _Nucleic Acids Research_ 32, pD493–D496, doi: https://doi.org/10.1093/nar/gkh103

[2] Pardiñas et al. (2018), "Common schizophrenia alleles are enriched in mutation-intolerant genes and in regions under strong background selection", _Nature Genetics_ 50, p381-389, doi: https://doi.org/10.1038/s41588-018-0059-2

### Glossary
* **SNP** - Single base-pair change that occurs on a specific genome position
* **rs ID** - Reference SNP ID, identfier for SNPs that map to identical genomic position, assigned by NCBI to reduce redundancy
