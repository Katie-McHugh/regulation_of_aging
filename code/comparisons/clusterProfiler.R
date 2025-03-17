###########################################################################
## GO-analysis using clusterProfiler ##
###########################################################################

#------------------------------------------------------------------------------
## Load in packages
#------------------------------------------------------------------------------

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# 
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)

# browseVignettes("clusterProfiler")
# BiocManager::install("topGO")
# browseVignettes("topGO")

BiocManager::install("org.Sc.sgd.db")

# Install biomaRt package (if not already installed)
BiocManager::install("biomaRt")

# Load the biomaRt library
library(biomaRt)
#------------------------------------------------------------------------------
## Load in gene lists #make sure using only lists of genic variants
#------------------------------------------------------------------------------

genes_dna<-read.csv("temp/genome/genic_GO_list.txt")
genes_rna<-read.csv("temp/transcriptome/DGEgene_list.txt")
genes_both<-read.csv("temp/comparisons/genic_combined_GO_list.txt")

#------------------------------------------------------------------------------
## Convert Gene Lists to EntreZ IDs
#------------------------------------------------------------------------------

# Connect to the Ensembl BioMart database
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "scerevisiae_gene_ensembl")

# You can check the available attributes (fields) and filters (options) in the Mart
attributes <- listAttributes(ensembl)
list(attributes)
filters <- listFilters(ensembl)
View(filters)

####### DNA list ###############
results_dna <- getBM(
  attributes = c("external_gene_name", "entrezgene_id"),  # Gene names and Entrez IDs
  filters = "external_gene_name",                      # Use gene names as filter
  values = genes_dna,                                 # Your list of gene names
  mart = ensembl                                       # The Ensembl mart
)

attributes <- listAttributes(ensembl)
head(attributes)
# Print the results
print(results_dna)

View(attributes)

dna_list<-results_dna$entrezgene_id

write.table(dna_list, file="temp/genome/dna_entrezIDs.txt", 
            row.names= FALSE, col.names= FALSE)

####### RNA list #####################

results_rna <- getBM(
  attributes = c("external_gene_name", "entrezgene_id"),  # Gene names and Entrez IDs
  filters = "external_gene_name",                      # Use gene names as filter
  values = genes_rna,                                 # Your list of gene names
  mart = ensembl                                       # The Ensembl mart
)

print(results_rna)

rna_list<-results_rna$entrezgene_id

write.table(rna_list, file="temp/transcriptome/rna_entrezIDs.txt", 
            row.names= FALSE, col.names= FALSE)

####### RNA list #########################

results_both <- getBM(
  attributes = c("external_gene_name", "entrezgene_id"),  # Gene names and Entrez IDs
  filters = "external_gene_name",                      # Use gene names as filter
  values = genes_both,                                 # Your list of gene names
  mart = ensembl                                       # The Ensembl mart
)

both_list<-results_both$entrezgene_id

write.table(both_list, file="temp/comparisons/rna&dna_entrezIDs.txt", 
            row.names= FALSE, col.names= FALSE)


#------------------------------------------------------------------------------
## Try to use clusterProfiler
#------------------------------------------------------------------------------

genes_dna<-as.list(genes_dna)
str(genes_dna)

?compareCluster
ck <- compareCluster(geneClusters = genes_dna, 
                     fun = enrichGO, 
                     OrgDb = "org.Sc.sgd.db")


ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 
#------------------------------------------------------------------------------
