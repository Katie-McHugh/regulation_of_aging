
##################################################################################
### Notes on Analysis comparing genomic and transcriptomic sequence data
##################################################################################

### Table of Contents: 
#---------------------------------------------------------------------------------------------------------------------------------------
1) 00_common_gene_lists.R
2) 00_local_variants.R
3) 02_Fig2_GO_terms.R
4) 01_genic_vs_nongenic.R
5) ) 01_SNP_annotation_comparisons.R
6) 02_FishersTest_annotations.R
7) 02_Figure4_barchart.R
8) 03_rna_supp_table.R

#---------------------------------------------------------------------------------------------------------------------------------------
### 1) 00_common_gene_lists.R
#---------------------------------------------------------------------------------------------------------------------------------------

#### Creates simplified lists of significant genes from genomic and transcriptomic analysis to use for comparisons

#---------------------------------------------------------------------------------------------------------------------------------------
### 2) 00_local_variants.R
#---------------------------------------------------------------------------------------------------------------------------------------

#### Script that is used to create Table 1, showing potential local cis-regulatory variants.  Can adjust the distance (default 10kb) used to determine what a "local" variant is.

#---------------------------------------------------------------------------------------------------------------------------------------
### 3) 02_Fig2_GO_terms.R
#---------------------------------------------------------------------------------------------------------------------------------------

#### Script that uses the outputs of the GO-term analyses performed using the Saccharomyces Genome Database GO analysis tool to visualize GO-term enrichment. Requires output of 01_genic_vs_nongenic.R to run GO-terms analysis.

#---------------------------------------------------------------------------------------------------------------------------------------
### 4) 01_genic_vs_nongenic.R
#---------------------------------------------------------------------------------------------------------------------------------------

#### Finds genomic regions that are implicated in both genomic and transcriptomic datasets

## generates list of genic terms for GO-term analysis

#---------------------------------------------------------------------------------------------------------------------------------------
### 5) 01_SNP_annotation_comparisons.R
#---------------------------------------------------------------------------------------------------------------------------------------
#### Script reorganized annotated SNP data to look for evidence of overrepresented categories in significant SNPs - output us used for Fishers Test in 02_Fisherstest_annotations.R

#---------------------------------------------------------------------------------------------------------------------------------------
### 6) 02_FishersTest_annotations.R
#---------------------------------------------------------------------------------------------------------------------------------------
#### Script that performs a Fisher's Exact test to find the probability of our observed proportion of annotation types in the significant SNP data given the overall proportions of each annotation type in the genome as a whole.  Dataset used is generated using 01_SNP_annotation_comparisons.R.
#---------------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------------------
### 7) 02_Figure4_barchart.R
#---------------------------------------------------------------------------------------------------------------------------------------
#### Script that creates Figure 4 annotation barchart
#---------------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------------------
### 8) 03_rna_supp_table.R
#---------------------------------------------------------------------------------------------------------------------------------------

#### Generates table showing relationships between DEGs and significant gene variants from genome data ("temp/comparisons/sig_FIRSTannotation.csv") and gene list from 00_common_gene_lists.R.  Preliminary version of suoplementary table 6

