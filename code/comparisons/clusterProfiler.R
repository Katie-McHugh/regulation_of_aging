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

# Install the necessary package
BiocManager::install("org.Sc.eg.db")
# Load the package
library(org.Sc.eg.db)





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

#-----------------------------------
####### DNA list ###############
#-----------------------------------
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

#-----------------------------------
####### RNA list ###################
#-----------------------------------

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

#-----------------------------------
####### Combined list ##############
#-----------------------------------

results_both <- getBM(
  attributes = c("external_gene_name", "entrezgene_id"),  
  filters = "external_gene_name",                      
  values = genes_both,                                 
  mart = ensembl                                      
)

both_list<-results_both$entrezgene_id

write.table(both_list, file="temp/comparisons/rna&dna_entrezIDs.txt", 
            row.names= FALSE, col.names= FALSE)


#------------------------------------------------------------------------------
## Try to use clusterProfiler ##not really working yet ##
#------------------------------------------------------------------------------

lists<-list(DNA=as.character(dna_list), 
            RNA=as.character(rna_list), 
            Combined=as.character(both_list))
            
str(lists)

# ck <- compareCluster(geneClusters = lists, 
#                      fun = enrichGO, 
#                      OrgDb = "org.Sc.sgd.db")

# ck <- setReadable(ck, OrgDb = "org.Sc.sgd.db", keyType="ENTREZID")


# dotplot(
#   ck,
#   x = "GeneRatio",
#   color = "p.adjust",
#   showCategory = 3,
#   size = NULL,
#   split = NULL,
#   font.size = 12,
#   title = "",
#   label_format = 30,
# )
dna_list<-as.character(dna_list)
head(dna_list)
ggo <- groupGO(gene     = dna_list,
               OrgDb    = "org.Sc.sgd.db",
               ont      = "CC",
               level    = 3,
               readable = TRUE)
View(ggo)


#------------------------------------------------------------------------------
# Try to load in GO-term data
#------------------------------------------------------------------------------
go_c<-read.table("results/GO_analyses/Genic/GO_genic_component.txt", 
                 header=TRUE, sep = "\t", quote = "")

go_p <- read.table("results/GO_analyses/Genic/GO_genic_process.txt", 
                header = TRUE, sep = "\t", quote = "")

#------------------------------------------------------------------------------
# Rename columns to match clusterProfiler
#------------------------------------------------------------------------------
temp_plot<-ggplot(go_c, aes(x = reorder(TERM, -CORRECTED_PVALUE), y = CORRECTED_PVALUE, size = NUM_LIST_ANNOTATIONS, color = CORRECTED_PVALUE)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "GO Terms", y = "P-value", title = "GO Term Enrichment") +
  theme_minimal() +
  coord_flip() + # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot(temp_plot)


