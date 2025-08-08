# regulation_of_aging


---------------------------------------------------------------------------------------------------------------------------------------

This repository contains the code and data for the analysis of: 

## Age-specific variation in genomic and transcriptomic signals reveal little evidence for cis-regulatory interactions modulating aging in Saccharomyces cerevisiae

### Repository content: 

#### [Data Folder](/data): 

This folder contains three subdirectories: [raw](/data/raw) , [clean](/data/clean) , and [design_files](/data/design_files).

##### [Raw](/data/raw)
1) filtered_snps.txt contains the "raw" SNP table that was generated using the Burke Lab pipeline
2) filtered_indels.txt contains the "raw" indel table that was generated using the Burke Lab pipeline
3) annotated_snps.txt contains the SNPeff annotations for the generated SNP table
4) annotated_indels.txt contains the SNPeff annotations for the generated indel table
5) founder_states.txt contains the haplotype frequencies for the data
6) gene_count_matrix.txt contains the gene count matrix from the transcriptomic sequencing (notes for generation are contained in the README.md file in the [transcriptomics code folder](/code/transcriptome).
7) all_founders.txt contains whole genome haplotype data for all founders (including some not used in this study)

##### [Clean](/data/clean)
1) scarsdata.txt contains the compiled bud scar count data 
2) rnaseq_results_batch_sigs01_edited.csv contains the results of the transcriptomic analysis (adjusted for the batch effect) and filtered for adjusted p<0.1.
3) GWAS_SNPS_cov20_maf05.csv contains the SNP table after filtering for a coverage of 20 and a minor allele frequency of 0.05
4) indels_cov20_maf05.csv contains the indel table after filtering for a coverage of 20 and a minor allele frequency of 0.05
5) CMHtest_cov_freq_matrix.csv
6) rnaseq_GO_terms.csv 
7) founder_states.txt gives the allele frequencies for each of the four isogenic founders and the 4S ancestor, and was used for haplotype analysis (generated from all_founders.txt using 00_founder_info.R script)


##### [Design Files](/data/design_files)
1) design.txt contains replicate identifying information for the transcriptomic analysis


#### [Code Folder](/code)
 This folder is broken up into three subfolders: 
 
1) [**genome**](/code/genome): Code for genomic variant analysis, haplotype analysis
 
2)  [**transcriptome**](/code/transcriptome): Code for differential gene expression analyses
 
3) [**comparisons**](/code/comparisons): Code to compare genomic and transcriptomic data and look into potential regulatory interactions
 
More information on the code within each folder is contained within the [code folder README](/code/README.md)

#### [Figures Folder](/figures)

This folder contains figures used in the manuscript. These include: 
1) _Fig1_Manhattan&Haps.pdf_: Figure 1: A manhattan plot of the nuclear genome showing differences in SNPs and indels between young and aged replicates
2) _Figure2_GO.jpeg_: Figure 2, a visualization of significant GO-terms
3) _Figure3_Volcano.pdf_ Figure 3: A Volcano plot of gene expression data, with significant gene variants labelled
4) _Figure4_annotations_barchart_color.pdf_ : Figure 4, which visualizes annotations from significant gene variants
**supplementary figure 1 is an image file of the sort**
5) _SuppFig2_BUDSCARS.odf_: Supplementary Figure 2, showing bud scar counts for each founder and the ancestor before and after sorting
6) _supp_fig_3_indels&SNPs_manhattan.pdf_ : Supplementary Figure 3, manhattan plot of mitochondrial indels and SNPs
7) _Supp_Fig4_Delta.pdf_ : Supplementary Figure 4, plot visualizing differences in haplotype frequencies across the genome
8) _SuppFig_heatmap_v2.pdf_: Supplementary Figure 5, heatmap visualizing DEGs


#### [Tables Folder](/tables)

This folder contains tables used in the manuscript generated in R. These include: 

1) _Table1_local_variants.csv_ : Table 1 in manuscript
2) _SuppTable2_GenomicVariantAnn.csv_ : Supplementary Table 2 in Manuscript, shows significant genetic variants from GWAS/CMH test
3) _ST3_GWAS_cov.txt_ : Supplementary Table 3, shows average coverage values for each replicate across the genome for SNPs and indels3) _ST6_variant_dist_table.csv_ : Supplementary Table 6, shows DEGs with additional info on their distance from significant gene variants (taken from SGD, SNPeff, and the UCSC genome database)

** Other supplementary tables were created directly in excel, using data from:  
- Supplementary Table 1:  budscars.R script in the [code folder](/code/budscars)
-  Supplementary Table 4 : GO-analysis results in the [results folder](/results/GO_analysis)
- Supplementary Table 5 : haplotype analysis 
- Supplementary Table 6 : transcriptomic analysis (see transcriptomics README)**

#### [Results Folder](/results)

This folder contains a subdirectory [GO_analyses](/results/GO_analyses) that contains multiple other subdirectories with the output from each gene ontology analyses, in addition to two directories for results associated with the genome data ([genome](/results/genome) ) and the transcriptome data ([transcriptome](/results/transcriptome) ). 

##### GO_analyses

###### Genic - 
1) GO_genic_component.txt contains the results for a GO-term analysis conducted on the significant SNPs that were located within genic regions in our GWAS analysis, looking at component with an FDR correction and 0.1 cutoff
2) GO_genic_process.txt contains the results for a GO-term analysis conducted on the significant SNPs that were located within genic regions in our GWAS analysis, looking at process with an FDR correction and 0.1 cutoff

**there is no GO_genic_function.txt because no functional components were implicated**


#### [Temp Folder](/temp)

This folder contains temporary files used for downstream analysis.  Each file will be generated using earlier code scripts, so these don't need to be downloaded.
 

