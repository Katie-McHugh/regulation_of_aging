#############################################################################
### Heatmap V3: Panels by function # complex heatmap
#############################################################################

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)
library(dplyr)

### Generating DGE Heatmap
#############################################################################

#### Note about LFC figure legend: 
#-------------------------------------------------------------------------------
## Pheatmap's side color legend (the LFC bar) is an estimate, and by default it
# will underestimate the labels assigned to the legend for the colors it uses
# e.g., the legend doesn't capture the full range of values 
# there is no way to manually correct pheatmap to tell it to re-assign 
# the min and max to the ACTUAL min and max--so the scale of the LFC bar 
# matches the data, but the actual values aren't perfect

## Spent a long time trying to append an accurate legend, but it doesn't seem to
# be possible using pheatmap
#-------------------------------------------------------------------------------
library(pheatmap)
## Load and Organize data

### Load normalized counts for visualization
norm_dds<-read.csv("temp/transcriptome/normalized_counts_deseq.csv")

### Load list of genes for headmap
selected_genes_adj<-read.csv("temp/transcriptome/RNA_genes_p<0.1.csv", row.names= "X")
selected_genes_adj2<-read.csv("temp/transcriptome/RNA_genes_p<0.05.csv", row.names= "X")

### load in design file
colData<-read.table("data/design_files/design.txt", header=TRUE, row.names = "sample")

# change "old" to "aged"
colData <- colData %>%
  mutate(across(where(is.character), ~ gsub("old", "aged", .)))

### Load in gene key
key<-read.table("temp/transcriptome/key_geneIDtoName.txt")

### subset norm_dds to significant genes
sel<-as.data.frame(selected_genes_adj$Gene_ID)
colnames(sel)[colnames(sel) == "selected_genes_adj$Gene_ID"] <- "Gene_ID"
norm_sigs<-merge(sel, norm_dds, by.x="Gene_ID", by.y="X", all.x=TRUE)

# Set key for old vs young replicates
age_colors<- c( "aged" = "black", "young"= "grey")
age_cols<-list(condition = age_colors)

### Convert gene ID to gene name, and make that the row name
norm_sigs2<-merge(norm_sigs, key, by="Gene_ID", all.x = TRUE)

## rename using SGD
norm_sigs2$Gene_Name[norm_sigs2$Gene_ID == "YLR050C"] <- "EMA19" ## Name recognized in SGD
norm_sigs2$Gene_Name[norm_sigs2$Gene_ID == "YCR015C"] <- "CTO1" ## Name recognized in SGD

rownames(norm_sigs2) <- norm_sigs2$Gene_Name

# Remove the non-numeric columns
norm_sigs2$Gene_ID <- NULL
norm_sigs2$Gene_Name <- NULL

# scale and convert to matrix
norm_adj<-t(scale(t(norm_sigs2))) 
#scale so each feature has the same mean/variance for visualization purposes 
# in the heatmap (prevents features with higher overall expression from washing 
# out signals of lower expression features)
norm_adj<-as.data.frame(norm_adj)
norm_adj<-unique(norm_adj)
norm_adj_mat<-as.matrix(norm_adj)
mode(norm_adj_mat) <- "numeric"
View(norm_adj_mat)

## Colors to distinguish between replicate pairs
subject_colors <- c(
  "pair1" = "#E41A1C", "pair2" = "#377EB8", "pair3" = "#4DAF4A", 
  "pair4" = "#984EA3", "pair5" = "#FF7F00", "pair6" = "#FFFF33", 
  "pair7" = "#A65628", "pair8" = "#F781BF", "pair9" = "#999999", 
  "pair10" = "#66C2A5", "pair11" = "#FC8D62", "pair12" = "#8DA0CB"
)

# Set up annotations for heat map
View(colData)
annotation_col <- data.frame(condition = colData$condition, pair=colData$subject)
rownames(annotation_col) <- rownames(colData) 
head(annotation_col) #match replicate label to age

### LOG FOLD CHANGE ANNOTATION
selected_genes_adj$logBM<-log10(selected_genes_adj$baseMean)
logfc_info<-unique(selected_genes_adj[,c("log2FoldChange", "logBM", "padj", "Gene_ID", "Gene_Name")])
rownames(logfc_info) <- logfc_info$Gene_Name

# Create a color palette for the LFC values
lfc_colors <- colorRampPalette(c("purple", "white", "darkgreen"))(100)
# bM_colors<- colorRampPalette(c("white", "black"))(100)



# Create a data frame for annotation, including LFC and BM
annotation_row <- data.frame(
  LFC = logfc_info$log2FoldChange 
  #, LBM = logfc_info$logBM
)

# extract just the logFC info
rownames(annotation_row) <- logfc_info$Gene_Name # Ensure row names match the heatmap data

head(annotation_row)
head(annotation_col)

# Define color scale for LFC values
annotation_colors <- list(
  condition = age_colors,
  pair= subject_colors,
  LFC = lfc_colors #,  # Use your color palette for LFC
  # LBM = bM_colors
)

### check that everything looks right
all(rownames(annotation_row) == rownames(norm_adj_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(norm_adj_mat))

View(logfc_info)

#-------------------------------------------------------------------------------
## Plot heatmap
### p < 0.1

heatmap<-pheatmap(
  norm_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_legend = TRUE, 
  fontsize_row = 10, # Adjust the font size of row names
  fontsize_col = 10, # Adjust the font size of column names
  cellheight = 12)
heatmap

#-------------------------------------------------------------------------------
## Convert to complex heatmap

c_heatmap<-ComplexHeatmap::pheatmap(
  norm_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_legend = TRUE, 
  fontsize_row = 10, 
  fontsize_col = 10,
  cellheight = 12)

ComplexHeatmap::pheatmap(
  norm_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_legend = TRUE, 
  fontsize_row = 10, 
  fontsize_col = 10,
  cellheight = 12)
print(c_heatmap)


pdf("test_heatmap.pdf", width=8, height=10)
draw(c_heatmap)
dev.off()
quartz()
draw(c_heatmap)


#---------------------------------------------------------------------------
## subset by panel 

## List genes belonging to each category ## some are listed more than once
group_list<-list(
"dna"=c("CLB6","UBC13", "CDC45", "RFA3", "MHR1", "MAG1", "ACT1", "POL12"),
"ion"=c("ATX1", "ATX2", "AHP1", "FIT2", "VMA10", "DNM1", "CTO1", "ISD11"),
"biosynth"=c("MET17", "EKI1", "CYB5", "ERG25", "IPT1", "TMA17"),
"protein"=c("CPR6", "SOM1", "EMA19", "SSA1", "SUE1", "TOS2","WSC4", "VID25", "SPG5", "MKK1","RPN6"),
"trans"=c("AIR2", "OPI1", "BAS1", "RRP36", "DIA4", "GCR1", "TRM5", "TOD6", "PCL5", "SUI2", "MRP1", "NDT80"),
"other"=c("YRO2", "PER33", "ARI1", "YDR543C", "YNR066C", "IRC18", "JIP4", "EFM4", "YJL218W", "COA2", "ARO10", "SPS4", "SOK1", "CWP1","CCW12","YLR363W-A")
)

gene_group <- unlist(group_list)
group_labels <- rep(names(group_list), lengths(group_list))
names(group_labels) <- gene_group
head(gene_group)


### subset by panel
dna<- norm_sigs2[norm_sigs2$Gene_Name %in% dna_list, ]
ion<- norm_sigs2[norm_sigs2$Gene_Name %in% ion_list, ]
biosynth<-norm_sigs2[norm_sigs2$Gene_Name %in% biosynth_list, ]
protein<-norm_sigs2[norm_sigs2$Gene_Name %in% protein_list, ]
trans<-norm_sigs2[norm_sigs2$Gene_Name %in% trans_list, ]
structural<-norm_sigs2[norm_sigs2$Gene_Name %in% structural_list, ]
other<-norm_sigs2[norm_sigs2$Gene_Name %in% other_list, ]

