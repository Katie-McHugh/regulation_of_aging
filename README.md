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

##### [Clean](/data/clean)
1) scarsdata.txt contains the compiled bud scar count data 
2) rnaseq_results_batch_sigs01_edited.csv contains the results of the transcriptomic analysis (adjusted for the batch effect) and filtered for adjusted p<0.1.
3) GWAS_SNPS_cov20_maf05.csv contains the SNP table after filtering for a coverage of 20 and a minor allele frequency of 0.05
4) indels_cov20_maf05.csv contains the indel table after filtering for a coverage of 20 and a minor allele frequency of 0.05
5) CMHtest_cov_freq_matrix.csv
6) rnaseq_GO_terms.csv 

##### [Design Files](/data/design_files)
1) design.txt contains replicate identifying information for the transcriptomic analysis


#### [Code Folder](/code)
 This folder is broken up into three subfolders: 
 
1) [**genome**](/code/genome): Code for genomic variant analysis, haplotype analysis
 
2)  [**transcriptome**](/code/transcriptome): Code for differential gene expression analyses
 
3) [**comparisons**](/code/comparisons): Code to compare genomic and transcriptomic data and look into potential regulatory interactions
 
More information on the code within each folder is contained within the [code folder README](/code/README.md)

#### [Figures Folder](/figures)

**need to clean up**

This folder contains figures used in the manuscript. These include: 
1) _Figure 1_: A manhattan plot of the nuclear genome showing differences in SNPs and indels between young and aged replicates
3) _Figure2

#### [Tables Folder](/tables)

This folder contains tables used in the manuscript. These include: 
1) 

#### [Results Folder](/results)

This folder contains a subdirectory [GO_analyses](/results/GO_analyses) that contains multiple other subdirectories with the output from each gene ontology analyses.

##### GO_analyses

###### Genic - 
1) GO_genic_component.txt contains the results for a GO-term analysis conducted on the significant SNPs that were located within genic regions in our GWAS analysis, looking at component with an FDR correction and 0.1 cutoff
2) GO_genic_process.txt contains the results for a GO-term analysis conducted on the significant SNPs that were located within genic regions in our GWAS analysis, looking at process with an FDR correction and 0.1 cutoff

**there is no GO_genic_function.txt because no functional components were implicated**


#### [Temp Folder](/temp)

This folder contains temporary files used for downstream analysis.  Each file should be generated using earlier code scripts, so these don't need to be downloaded.
 
 #### [Temp Figures Folder](/temp_figures)
 
 This folder contains figures used for preliminary visualization
 
  #### [Temp Tables Folder](/temp_tables)
  
  This folder contains tables used for preliminary investigations
  
  #### [Manuscript Folder](/manuscript)
  
  This folder contains...
 

