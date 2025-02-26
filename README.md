# regulation_of_aging

## Genome- and transcriptome- wide association study of aging in S. cerevisiae. PhD dissertation chapter.
---------------------------------------------------------------------------------------------------------------------------------------

This repository contains the code and data for the analysis of: 

xxx

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
1) rnaseq_results_batch_sigs01_edited.csv

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
1) _Figure 1_: A manhattan plot of the nuclear genome showing differences in SNPs and indels between young and aged replicates

#### [Tables Folder](/tables)

This folder contains tables used in the manuscript. These include: 
1) 

#### [Results Folder](/results)

This folder contains a subfolder [GO_analyses](/results/GO_analyses) that contains the output from the gene ontology analyses.
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
 

