###########################################################################
## GO-analysis Visualization ##
###########################################################################

#------------------------------------------------------------------------------
## Load in packages
#------------------------------------------------------------------------------

 if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

 BiocManager::install("pathview")
 BiocManager::install("enrichplot")
 #library(enrichplot)

 #BiocManager::install("org.Sc.sgd.db")

 # Install biomaRt package (if not already installed)
 BiocManager::install("biomaRt")

 # Load the biomaRt library
 library(biomaRt)

 # Install the necessary package
 BiocManager::install("org.Sc.eg.db")
 # Load the package
 library(org.Sc.eg.db)


# #------------------------------------------------------------------------------
# ## Load in gene lists #make sure using only lists of genic variants
# #------------------------------------------------------------------------------
# 
# genes_dna<-read.csv("temp/genome/genic_GO_list.txt")
# genes_rna<-read.csv("temp/transcriptome/DGEgene_list.txt")
# genes_both<-read.csv("temp/comparisons/genic_combined_GO_list.txt")
# 
# #------------------------------------------------------------------------------
# ## Convert Gene Lists to EntreZ IDs
# #------------------------------------------------------------------------------
# 
# # Connect to the Ensembl BioMart database
# ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "scerevisiae_gene_ensembl")
# 
# # You can check the available attributes (fields) and filters (options) in the Mart
# attributes <- listAttributes(ensembl)
# filters <- listFilters(ensembl)
# 
# #-----------------------------------
# ####### DNA list ###############
# #-----------------------------------
# results_dna <- getBM(
#   attributes = c("external_gene_name", "entrezgene_id"),  # Gene names and Entrez IDs
#   filters = "external_gene_name",                      # Use gene names as filter
#   values = genes_dna,                                 # Your list of gene names
#   mart = ensembl                                       # The Ensembl mart
# )
# 
# attributes <- listAttributes(ensembl)
# head(attributes)
# # Print the results
# print(results_dna)
# 
# View(attributes)
# 
# dna_list<-results_dna$entrezgene_id
# 
# write.table(dna_list, file="temp/genome/dna_entrezIDs.txt", 
#             row.names= FALSE, col.names= FALSE)
# 
# #-----------------------------------
# ####### RNA list ###################
# #-----------------------------------
# 
# results_rna <- getBM(
#   attributes = c("external_gene_name", "entrezgene_id"),  # Gene names and Entrez IDs
#   filters = "external_gene_name",                      # Use gene names as filter
#   values = genes_rna,                                 # Your list of gene names
#   mart = ensembl                                       # The Ensembl mart
# )
# 
# print(results_rna)
# 
# rna_list<-results_rna$entrezgene_id
# 
# write.table(rna_list, file="temp/transcriptome/rna_entrezIDs.txt", 
#             row.names= FALSE, col.names= FALSE)
# 
# #-----------------------------------
# ####### Combined list ##############
# #-----------------------------------
# 
# results_both <- getBM(
#   attributes = c("external_gene_name", "entrezgene_id"),  
#   filters = "external_gene_name",                      
#   values = genes_both,                                 
#   mart = ensembl                                      
# )
# 
# both_list<-results_both$entrezgene_id
# 
# write.table(both_list, file="temp/comparisons/rna&dna_entrezIDs.txt", 
#             row.names= FALSE, col.names= FALSE)
# 

#------------------------------------------------------------------------------
# Try to load in GO-term data
#------------------------------------------------------------------------------
go_c<-read.table("results/GO_analyses/Genic/GO_genic_component.txt", 
                 header=TRUE, sep = "\t", quote = "")

go_p <- read.table("results/GO_analyses/Genic/GO_genic_process.txt", 
                header = TRUE, sep = "\t", quote = "")

go_both_p<-read.table("results/GO_analyses/Genic/GO_combined_process.txt", 
                 header=TRUE, sep = "\t", quote = "")

go_both_f <- read.table("results/GO_analyses/Genic/GO_function_combined.txt", 
                   header = TRUE, sep = "\t", quote = "")
go_both_c<- read.table("results/GO_analyses/Genic/genic_combined_component_GO.txt", 
                      header = TRUE, sep = "\t", quote = "")

## add column identifiers to combine lists
go_c$gene_list<-"DNA"
go_p$gene_list<-"DNA"
go_c$ann<-"component"
go_p$ann<-"process"
go_both_p$gene_list<-"both"
go_both_f$gene_list<-"both"
go_both_c$gene_list<-"both"
go_both_p$ann<-"process"
go_both_f$ann<-"function"
go_both_c$ann<-"component"




go_dna <- bind_rows(go_c, go_p, go_both_p, go_both_c, go_both_f)
View(go_dna)

#------------------------------------------------------------------------------
# Create GO-plots
#------------------------------------------------------------------------------

DNA_c_plot<-ggplot(go_c, aes(x = reorder(TERM, -CORRECTED_PVALUE), 
                             y = CORRECTED_PVALUE, size = NUM_LIST_ANNOTATIONS, 
                             color = CORRECTED_PVALUE)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "GO Terms", y = "P-value", title = "GO_dna_component") +
  theme_minimal() +
  coord_flip() + # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DNA_p_plot<-ggplot(go_p, aes(x = reorder(TERM, -CORRECTED_PVALUE), 
                             y = CORRECTED_PVALUE, size = NUM_LIST_ANNOTATIONS, 
                             color = CORRECTED_PVALUE)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "GO Terms", y = "P-value", title = "GO DNA process") +
  theme_minimal() +
  coord_flip() + # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot(DNA_p_plot)
plot(DNA_c_plot)

## combined plot

# Create the plot
c_plot<-ggplot(go_dna, aes(x = reorder(TERM, -CORRECTED_PVALUE), y = CORRECTED_PVALUE, 
                          size = NUM_LIST_ANNOTATIONS, color = CORRECTED_PVALUE)) +
  geom_point(alpha = 0.5) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "GO Terms", y = "Corrected P-value", 
       title = "GO Term Enrichment for Gene List 1 (Component & Function)") +
  theme_minimal() +
  coord_flip() +  # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ ann, scales = "free_y")  # Separate Component vs Function

plot(c_plot)

## same plot, different colors

# Create the plot
plot2<-ggplot(go_dna, aes(x = reorder(TERM, -CORRECTED_PVALUE), y = CORRECTED_PVALUE, 
                          size = NUM_LIST_ANNOTATIONS, color = gene_list)) +
  geom_point(alpha = 0.5) +  
  scale_color_manual(values = c("DNA" = "blue", "both" = "red")) +              # Custom colors for Component and Function
  labs(x = "GO Terms", y = "Corrected P-value", 
       title = "GO Term Enrichment") +
  theme_minimal() +
  coord_flip() +  # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(name = "Number of Annotated Genes")  +
  facet_wrap(~ ann, scales = "free_y")  
plot(plot2)
ggsave(filename = "figures/FigureXX_GO_sepType.jpeg", plot = plot2 , 
       width = 12, height = 4)


### different style

# relevel lists: 
go_dna$gene_list<-as.factor(go_dna$gene_list)
go_dna <- go_dna %>%
  mutate(gene_list = recode(gene_list, "both" = "Combined"))
go_dna$gene_list<-relevel(go_dna$gene_list, ref="DNA")

plot3<-ggplot(go_dna, aes(x = reorder(TERM, -CORRECTED_PVALUE), y = CORRECTED_PVALUE, 
                          size = NUM_LIST_ANNOTATIONS, color = ann)) +
  geom_point(alpha = 1) +  
  scale_color_manual(name = "GO-Term \nCategory", 
                     values = c("component" = "blue", 
                                "process" = "red",
                                "function"="goldenrod4")) +                      # Custom colors for Component and Function
  labs(x = "GO Terms", y = "Corrected P-value", 
       title = "GO Term Enrichment") +
  theme_minimal() +
  coord_flip() +  # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +              # Wrap y-axis labels
  scale_size_continuous(name = "Number of \nAnnotated \nGenes")  +
  facet_wrap(~gene_list, scales = "fixed")  + 
  theme_light() +
  theme(strip.text = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10))

plot(plot3)

library(RColorBrewer)
plot4 <- plot3 + scale_color_manual(values = brewer.pal(3, "Dark2"))
plot(plot4)

ggsave(filename = "figures/FigureXX_GO_sepList.jpeg", plot = plot4 , width = 8, height = 7)

### maybe try separating by process vs component vs etc with different gene 
## lists as the different colors
#------------------------------------------------------------------------------


### Version 2
###############################################################################
#### Use genes implicated in both lists instead of either list
###############################################################################

library(dplyr)
library(ggplot2)
library(stringr)

#------------------------------------------------------------------------------
# Try to load in GO-term data
#------------------------------------------------------------------------------
go_c<-read.table("results/GO_analyses/Genic/GO_genic_component.txt", 
                 header=TRUE, sep = "\t", quote = "")

go_p <- read.table("results/GO_analyses/Genic/GO_genic_process.txt", 
                   header = TRUE, sep = "\t", quote = "")

go_overlap_p<-read.table("results/GO_analyses/Overlap_List/Overlap_list_process.txt", 
                      header=TRUE, sep = "\t", quote = "")



## add column identifiers to combine lists
go_c$gene_list<-"DNA"
go_p$gene_list<-"DNA"
go_c$ann<-"component"
go_p$ann<-"process"
go_overlap_p$gene_list<- "both"
go_overlap_p$ann<-"process"


go_dna <- bind_rows(go_c, go_p, go_overlap_p)
View(go_dna)

#------------------------------------------------------------------------------
# Create GO-plots
#------------------------------------------------------------------------------

DNA_c_plot<-ggplot(go_c, aes(x = reorder(TERM, -CORRECTED_PVALUE), 
                             y = CORRECTED_PVALUE, size = NUM_LIST_ANNOTATIONS, 
                             color = CORRECTED_PVALUE)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "GO Terms", y = "P-value", title = "GO_dna_component") +
  theme_minimal() +
  coord_flip() +                                                                # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DNA_p_plot<-ggplot(go_p, aes(x = reorder(TERM, -CORRECTED_PVALUE),
                             y = CORRECTED_PVALUE, size = NUM_LIST_ANNOTATIONS, 
                             color = CORRECTED_PVALUE)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "GO Terms", y = "P-value", title = "GO DNA process") +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot(DNA_p_plot)
plot(DNA_c_plot)

## combined plot

# Create the plot
c_plot<-ggplot(go_dna, aes(x = reorder(TERM, -CORRECTED_PVALUE), y = CORRECTED_PVALUE, 
                           size = NUM_LIST_ANNOTATIONS, color = CORRECTED_PVALUE)) +
  geom_point(alpha = 0.5) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "GO Terms", y = "Corrected P-value", 
       title = "GO Term Enrichment for Gene List 1 (Component & Function)") +
  theme_minimal() +
  coord_flip() +                                                                # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ ann, scales = "free_y") +                                        # Separate Component vs Function
  theme_light() +
  theme(strip.text = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10))

plot(c_plot)

## same plot, different colors

#if one isn't showing up, the name is probably wrong...this is how to fix it
go_dna <- go_dna %>%
  mutate(gene_list = recode(gene_list, "both" = "Shared"))                       #if one isn'y showing up, the name is probably wrong...this is how to fix it

# Create the plot
plot2<-ggplot(go_dna, aes(x = reorder(TERM, -CORRECTED_PVALUE),
                          y = CORRECTED_PVALUE, 
                          size = NUM_LIST_ANNOTATIONS, color = gene_list)) +
  geom_point(alpha = 0.5) +  
  scale_color_manual(values = c("DNA" = "blue", "Shared" = "red")) +             # Custom colors for Component and Function 
  labs(x = "GO Terms", y = "Corrected P-value", 
       title = "GO Term Enrichment") +
  theme_minimal() +
  coord_flip() +  # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(name = "Number of \nAnnotated Genes")  +
  facet_wrap(~ ann, scales = "free_y", dir = "v")  +                             #dir v stacks vertically instead of horizontally
  theme_light() +
  theme(strip.text = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10))

plot(plot2)
ggsave(filename = "figures/FigureXX_GO_sepType_v3.jpeg", plot = plot2 ,
       width = 7, height = 7) 


### different style

# relevel lists: 
go_dna$gene_list<-as.factor(go_dna$gene_list)
go_dna <- go_dna %>%
  mutate(gene_list = recode(gene_list, "both" = "Shared"))
go_dna$gene_list<-relevel(go_dna$gene_list, ref="DNA")

plot3<-ggplot(go_dna, aes(x = reorder(TERM, -CORRECTED_PVALUE), y = CORRECTED_PVALUE, 
                          size = NUM_LIST_ANNOTATIONS, color = ann)) +
  geom_point(alpha = 1) +  
  scale_color_manual(name = "GO-Term \nCategory", 
                     values = c("component" = "blue", 
                                "process" = "red", 
                                "function"="goldenrod4")) +                      # Custom colors for Component and Function
  labs(x = "GO Terms", y = "Corrected P-value", 
       title = "GO Term Enrichment") +
  theme_minimal() +
  coord_flip() +  # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +               # Wrap y-axis labels
  scale_size_continuous(name = "Number of \nAnnotated \nGenes")  +
  facet_wrap(~gene_list, scales = "fixed")  + 
  theme_light() +
  theme(strip.text = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10))

plot(plot3)

library(RColorBrewer)
plot4 <- plot3 + scale_color_manual(values = brewer.pal(3, "Dark2"))
plot(plot4)

ggsave(filename = "figures/FigureXX_GO_sepList_v2.jpeg", plot = plot4 , width = 8, height = 7)

### maybe try separating by process vs component vs etc with different gene 
## lists as the different colors

