# README for Code Folder

## Notes on code folder contents: 

### Intro files: 

1) README (this file)
2) packages.R contains all the libraries and packages needed to run scripts in this folder (including the following sub-folders)

### Sub-folders: 

1) [genome](/code/genome) contains all scripts needed for genome-wide analysis
2) [transcriptome](/code/transcriptome) folder contains all scripts needed for transcriptome analysis
3) [comparisons](/comparisons) folder contains scripts for regulatory analysis and other combined analyses

#### [genome](/code/genome) folder: 

1) SNPs_filter_maf&cov.R contains the script used to filter the preliminary SNP table for minor allele frequencies and coverage. 
2) indels_filter_maf&cov.R contains the script used to filter the preliminary indel table for minor allele frequencies and coverage.
3) CMHtest_SNPdata.R contains the script used to run the CMH test for the SNP data
4) indels_cmhtest.R contains the script used to run the CMH test for the indel data
5) nuclearSNPs_ManhattanPlot.R is the script used to generate a manhattan plot for the nuclear SNPs
6) mitoSNPs_ManhattanPlot.R is the script used to generate a manhattan plot for the mitochondrial SNPs
7) indels_nuclear_manhattan.R is the script used to generate a manhattan plot for the nuclear indels
8) indels_mito_manhattan.R is the script used to generate a manhattan plot for the mitochondrial indels
9) indels&snps_nuclear_manhattan.R is the script used to generate a manhattan plot that overlays the indels on top of the SNP manhattan plot - incomplete/needs to be fixed.
10) Haplotype_estimation_eNotebook_version.R is the script that was used to estimate the haploype frequencies for each replicate. This script was transferred over after use, and was used iteratively to run for each chromosome and then they were all bound together into a single file. 
11) haplotype_differences.R is the script used to calculate the difference in haplotype frequencies across the genome to use for plotting
12) supp_table_sig_variant_list.R is the script used to create the supplementary table of all significant variants from teh genome data

#### [transcriptome](/code/transcriptome) folder:

1) batch_adjust_DESEQ.R is the script that takes the gene count matrix generated in the description below and performs a batch correction based on the sorting date, performs some additional filtering to remove genes with low mapping, and creates a DESEQ object for later analysis
2) DGE_gene_lists.R is the script used to create the significant gene lists for GO-term analysis, and to be used to compare against the genome data
3) heatmap.R is the script used to create the heatmap figures for the transcriptome data
4) rna_supp_table.R is the script used to create the supplementary data of significant transcripts and their distances from implicated genome variants
5) README.md file that contains notes on data filtering/cleaning prior to analysis

#### [comparisons](/comparisons) folder: 
1) 00_common_gene_lists.R - this script should be run before running any other scripts in this folder. It takes the significant gene lists from the genomic data and extracts the appropriate annotations for use in later analysis.
2) dna_rna_gene_comparison.R uses files generated from 00_common_gene_lists.R to find genes that are implicated in both genomic and transcriptomic analyses
3) genic_vs_nongenic.R uses files generated from 00_common_gene_lists.R to compare which variants are found within coding vs non-coding regions, and to create a combined list of genic significant variants for genomic and transcriptomic GO-analysis.
4) annotation_comparisons.R uses files generated from 00_common_gene_lists.R to compare the frequency of annotation categories observed in our significant genomic data to the entire annotated genome.
5) annotations_pie_chart.R uses files generated from annotation_comparisons.R to generate a pie chart that visualized the annotations (Figure)
6) FishersTest_annotations.R performs a statistical test on the annotation proportions from annotation_comparisons.R


... to be continued...


#### outputs should be dumped into temp folders

