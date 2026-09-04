#  Notes on code folder contents: 

## Intro files: 

1) README (this file)
2) packages.R contains all the libraries and packages needed to run scripts in this folder (including the following sub-folders)

##  Sub-folders: 

1) [genome](/code/genome) contains all scripts needed for whole-genome sequencing analysis
2) [transcriptome](/code/transcriptome) folder contains all scripts needed for transcriptome analysis
3) [comparisons](/comparisons) folder contains scripts for regulatory analysis and other combined analyses
4)  [**budscars**](/code/budscars) folder contains script for bud scar phenotyping analysis

### [genome](/code/genome) folder: 

**number prefixes indicate sequential scripts**

1) 00_SNPs_filter_maf&cov.R contains the script used to filter the preliminary SNP table for minor allele frequencies and coverage. Output sent to "clean" folder as GWAS_SNPS_cov20_maf05.csv.
2) 00_indels_filter_maf&cov.R contains the script used to filter the preliminary indel table for minor allele frequencies and coverage. Output sent to "clean" folder as indels_cov20_maf05.csv.
3) 00_founder_info.R contains code to generate a table of information for the founder haplotypes using founder_states.txt file (generated from an earlier experiment). 
4) 01_CMHtest_SNPdata.R contains the script used to run the CMH test on GWAS_SNPS_cov20_maf05.csv for the SNP data.  Output is sent to results folder.
5) 01_indels_cmhtest.R contains the script used to run the CMH test for the indel data, output to results folder.
6) 01_Haplotype_estimation_eNotebook_version.R is the script that was used to estimate the haploype frequencies for each replicate. This script was transferred over to repository after use. It was used iteratively to run for each chromosome and then they were all bound together into a single file for each founder (founderXX_10kb_cov20_maf5.txt) and deposited in the "temp" folder.
7) 01_freq&cov_matrix.R contains code to generate a table of frequency and coverage values at each SNP position used to create plots in the 02_freq&cov_plots.R script. It also creates supplementary Table 2 of GWAS coverage averages. 
8) 02_Delta_plots.R contains the script used to generate supplementary figure 8 (measures of haplotype differentiation)
9) 02_SuppFig3_indels&snps_mito_manhattan.R is the script used to generate supplementary figure 3, a manhattan plot that overlays the indels on top of the SNP manhattan plot for the mitochondria
10) 02_supp_table_sig_variant_list.R is the script used to create the supplementary table 2 containing all significant variants from the genome data (a small edit was later done to fill in annotations for the 2 indels that are not annotated in SNPeff)
11) 02_Figure1.R is the script used to generate figure 2, containing the nuclear SNP and indel manhattan plot and founder haplotype frequencies across the genome.
12) SFS.R is a script to generate site frequency spectra for Supplementary Figure 5
13) Coverage_figures.R is a script to generate plots to visualize coverage, including Supplementary Figure 6. 

#### [Extras](/code/genome/Extras) folder: 

**contains code not used in final manuscript, but previously used in other presentation formats** 

1) 02_nuclearSNPs_ManhattanPlot.R is the script used to generate a manhattan plot for the nuclear SNPs
2) 02_haplotype_differences.R is the script used to calculate the difference in haplotype frequencies across the genome to use for plotting
3) 02_mitoSNPs_ManhattanPlot.R is the script used to generate a manhattan plot for the mitochondrial SNPs
4) 02_indels_nuclear_manhattan.R is the script used to generate a manhattan plot for the nuclear indels
5) 02_indels_mito_manhattan.R is the script used to generate a manhattan plot for the mitochondrial indels
6) 02_indels&snps_nuclear_manhattan.R is the script used to generate a manhattan plot that overlays the indels on top of the SNP manhattan plot
7) 02_freq&cov_plots.R script used to plot coverage and frequency by replicate for individual SNPs.  Uses table generated in 01_freq&cov_matrix.R.


### [transcriptome](/code/transcriptome) folder:

1) 01_batch_adjust_DESEQ.R is the script that takes the gene count matrix generated in the description below and performs a batch correction based on the sorting date, performs some additional filtering to remove genes with low mapping, and creates a DESEQ object for later analysis. Also creates supplementary figure 3B&D (PCA).
2) 02_DGE_gene_lists.R is the script used to create the significant gene lists for GO-term analysis, and to be used to compare against the genome data
3) 03_heatmap.R is the script used to create the heatmap figures for the transcriptome data
4) 03_uncorrected_heatmap.R is a script visualizing the same heatmap, but without ComBatseq batch correction (used in supplementary figure)
 5) 03_Volcano_plot.R creates a volcano plot to show p-values and log2FC between DEGs using adjusted counts matrix from  01_batch_adjust_DESEQ.R and and gene list and key from 02_DGE_gene_lists.R.
6) README.md file that contains notes on data filtering/cleaning prior to analysis

Extras: 
Unused folder contains test scripts that were unfinished/scrapped (breaking up heatmap into panels)

### [comparisons](/comparisons) folder: 
1) 00_common_gene_lists.R - this script should be run before running any other scripts in this folder. It takes the significant gene lists from the genomic data and extracts the appropriate annotations for use in later analysis.
2) 00_local_variants.R use a table (results/variant_dist_table.csv) generated by combining the RNA-seq results with gene position data collected from the Saccharomyces genome database to find local variants near DEGs
3) 01_genic_vs_nongenic.R uses files generated from 00_common_gene_lists.R to compare which variants are found within coding vs non-coding regions, and to create a combined list of genic significant variants for genomic and transcriptomic GO-analysis.
4) 01_SNP_annotation_comparisons.R reorganized annotated SNP data to look for evidence of overrepresented categories in significant SNPs - output us used for Fishers Test in 02_Fisherstest_annotations.R
5) 02_Fig2_GO_terms.R uses results of GO-term analysis (which requires output from 01_genic_vs_nongenic.R) to create Figure 2 to visualize GO terms.
6) 02_FishersTest_annotations.R performs a statistical test on the annotation proportions from 01_annotation_comparisons.R
7) 02_Figure4_barchart.R creates Figure 4 annotation barchart
8) 03_rna_supp_table.R generates table showing relationships between DEGs and significant gene variants from genome data ("temp/comparisons/sig_FIRSTannotation.csv") and gene list from 00_common_gene_lists.R.  Preliminary version of supplementary table 6

#### [Extras](/code/comparisons/Extras) folder: 
9) 01_dna_rna_gene_comparison.R uses files generated from 00_common_gene_lists.R to find genes that are implicated in both genomic and transcriptomic analyses
10) 01_annotation_comparisons.R uses files generated from 00_common_gene_lists.R to compare the frequency of annotation categories observed in our significant genomic data to the entire annotated genome.
11) 00_local_variants_1kb.R same as local_variants.R, just with a more restrictive (1kb) range
12) coding_vs_noncoding_fig.R creates an infographic representing the number of coding variants, but does not require any of the previous scripts or data.
13) 02_annotations_pie_chart.R uses files generated from 01_annotation_comparisons.R to generate a pie chart that visualized the annotations (Figure)
14) 02_annotation_comparisons_barchart.R uses files generated from 01_annotation_comparisons.R create a bar chart figure that visualizes the annotations
15) 02_GO_term_plot.R uses the output of the GO-term analysis (which used the gene lists from 01_genic_vs_nongenic.R) to generate a GO-term figure representing overrepresented aspects
16) 02_DNA_RNA_distances.R creates a figure to show the number of differentially expressed genes with nearby significant gene variants.



### [budscars](/budscars) folder: 
1) budscars.R contains the script used to perform the statistical tests for budscar counts and to create Supplementary Table 1. This script also creates Supplementary Figure 2.

_outputs should be dumped into temp folders_

