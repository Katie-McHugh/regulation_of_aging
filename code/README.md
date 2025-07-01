#  Notes on code folder contents: 

## Intro files: 

1) README (this file)
2) packages.R contains all the libraries and packages needed to run scripts in this folder (including the following sub-folders)

##  Sub-folders: 

1) [genome](/code/genome) contains all scripts needed for genome-wide analysis
2) [transcriptome](/code/transcriptome) folder contains all scripts needed for transcriptome analysis
3) [comparisons](/comparisons) folder contains scripts for regulatory analysis and other combined analyses

### [genome](/code/genome) folder: 

**number prefixes indicate sequential scripts**

1) 00_SNPs_filter_maf&cov.R contains the script used to filter the preliminary SNP table for minor allele frequencies and coverage. 
2) 00_indels_filter_maf&cov.R contains the script used to filter the preliminary indel table for minor allele frequencies and coverage.
3) 00_founder_info.R contains code to generate a table of information for the founder haplotypes. 
4) 01_CMHtest_SNPdata.R contains the script used to run the CMH test for the SNP data
5) 01_indels_cmhtest.R contains the script used to run the CMH test for the indel data
6) 01_freq&cov_matrix.R contains code to generate a table of frequency and coverage values at each SNP position used to create plots in the 02_freq&cov_plots.R script.
7) 01_Haplotype_estimation_eNotebook_version.R is the script that was used to estimate the haploype frequencies for each replicate. This script was transferred over after use, and was used iteratively to run for each chromosome and then they were all bound together into a single file. 
8) 02_nuclearSNPs_ManhattanPlot.R is the script used to generate a manhattan plot for the nuclear SNPs
9) 02_mitoSNPs_ManhattanPlot.R is the script used to generate a manhattan plot for the mitochondrial SNPs
10) 02_indels_nuclear_manhattan.R is the script used to generate a manhattan plot for the nuclear indels
11) 02_indels_mito_manhattan.R is the script used to generate a manhattan plot for the mitochondrial indels
12) 02_indels&snps_nuclear_manhattan.R is the script used to generate a manhattan plot that overlays the indels on top of the SNP manhattan plot
13) 02_indels&snps_mito_manhattan.R is the script used to generate a manhattan plot that overlays the indels on top of the SNP manhattan plot for the mitochondria
14) 02_haplotype_differences.R is the script used to calculate the difference in haplotype frequencies across the genome to use for plotting
15) 02_supp_table_sig_variant_list.R is the script used to create the supplementary table of all significant variants from teh genome data
16) 02_freq&cov_plots.R script used to plot coverage and frequency by replicate for individual SNPs.  Uses table generated in 01_freq&cov_matrix.R.

### [transcriptome](/code/transcriptome) folder:

1) batch_adjust_DESEQ.R is the script that takes the gene count matrix generated in the description below and performs a batch correction based on the sorting date, performs some additional filtering to remove genes with low mapping, and creates a DESEQ object for later analysis
2) DGE_gene_lists.R is the script used to create the significant gene lists for GO-term analysis, and to be used to compare against the genome data
3) heatmap.R is the script used to create the heatmap figures for the transcriptome data
4) rna_supp_table.R is the script used to create the supplementary data of significant transcripts and their distances from implicated genome variants
 5) 03_Volcano_plot.R creates a volcano plot to show p-values and log2FC between DEGs using adjusted counts matrix from  01_batch_adjust_DESEQ.R and and gene list and key from 02_DGE_gene_lists.R.
6) README.md file that contains notes on data filtering/cleaning prior to analysis
7) Unused folder contains test scripts that were unfinished/scrapped (breaking up heatmap into panels)

### [comparisons](/comparisons) folder: 
1) 00_common_gene_lists.R - this script should be run before running any other scripts in this folder. It takes the significant gene lists from the genomic data and extracts the appropriate annotations for use in later analysis.
2) dna_rna_gene_comparison.R uses files generated from 00_common_gene_lists.R to find genes that are implicated in both genomic and transcriptomic analyses
3) genic_vs_nongenic.R uses files generated from 00_common_gene_lists.R to compare which variants are found within coding vs non-coding regions, and to create a combined list of genic significant variants for genomic and transcriptomic GO-analysis.
4) annotation_comparisons.R uses files generated from 00_common_gene_lists.R to compare the frequency of annotation categories observed in our significant genomic data to the entire annotated genome.
5) annotations_pie_chart.R uses files generated from annotation_comparisons.R to generate a pie chart that visualized the annotations (Figure)
6) FishersTest_annotations.R performs a statistical test on the annotation proportions from annotation_comparisons.R


... to be continued...


_outputs should be dumped into temp folders_

