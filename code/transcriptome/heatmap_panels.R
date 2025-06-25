#############################################################################
### Heatmap V2: Panels by function
#############################################################################

library(dplyr)
library(pheatmap)

## Load gene expression data
#---------------------------------------------------------------------------
### Load normalized counts for visualization
norm_dds<-read.csv("temp/transcriptome/normalized_counts_deseq.csv")

### Load list of genes for heatmap
selected_genes_adj<-read.csv("temp/transcriptome/RNA_genes_p<0.1.csv", row.names= "X")

### load in design file
colData<-read.table("data/design_files/design.txt", header=TRUE, row.names = "sample")

# change "old" to "aged"
colData <- colData %>%
  mutate(across(where(is.character), ~ gsub("old", "aged", .)))

### Load in gene key
key<-read.table("temp/transcriptome/key_geneIDtoName.txt")


#---------------------------------------------------------------------------
## Reformat for heatmap

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

#---------------------------------------------------------------------------
## subset by panel 

## List genes belonging to each category ## some are listed more than once
dna_repair<-c("CLB6","UBC13", "CDC45", "RFA3", "MHR1", "MAG1", "ACT1", "POL12")
ion<-c("ATX1", "ATX2", "AHP1", "FIT2", "VMA10", "DNM1", "CTO1", "ISD11")
biosynth<-c("MET17", "EKI1", "CYB5", "ERG25", "IPT1", "TMA17")
protein<-c("CPR6", "SOM1", "EMA19", "SSA1", "SUE1", "WSC4", "TMA17", "VID25", "SPG5", "RPN6")
trans_reg<-c("AIR2", "OPI1", "BAS1", "RRP36", "DIA4", "GCR1", "TRM5", "TOD6", "PCL5", "SUI2", "MRP1", "NDT80")
structure<-c("CWP1", "TOS2", "CCW12", "MKK1", "ACT1")
other<-c("YRO2", "PER33", "ARI1", "YDR543C", "YNR066C", "IRC18", "JIP4", "EFM4", "YJL218W", "COA2", "ARO10", "SPS4", "SOK1", "YLR363W-A")

norm_sigs2$Gene_ID <- NULL

### subset by panel
dna_rep<- norm_sigs2[norm_sigs2$Gene_Name %in% dna_repair, ]
ion_reg<- norm_sigs2[norm_sigs2$Gene_Name %in% ion, ]
biosynth2<-norm_sigs2[norm_sigs2$Gene_Name %in% biosynth, ]
protein_reg<-norm_sigs2[norm_sigs2$Gene_Name %in% protein, ]
trans_reg2<-norm_sigs2[norm_sigs2$Gene_Name %in% trans_reg, ]
structural<-norm_sigs2[norm_sigs2$Gene_Name %in% structure, ]
other2<-norm_sigs2[norm_sigs2$Gene_Name %in% other, ]

#---------------------------------------------------------------------------
# Remove the non-numeric columns
dna_rep$Gene_Name <- NULL
ion_reg$Gene_Name <- NULL
biosynth2$Gene_Name <- NULL
protein_reg$Gene_Name <- NULL
trans_reg2$Gene_Name <- NULL
structural$Gene_Name <- NULL
other2$Gene_Name <- NULL

#---------------------------------------------------------------------------
# DNA repair, replication, cell cycle regulation

# scale and convert to matrix
dna_adj<-t(scale(t(dna_rep))) 
#scale so each feature has the same mean/variance for visualization purposes 
# in the heatmap (prevents features with higher overall expression from washing 
# out signals of lower expression features)
dna_adj<-as.data.frame(dna_adj)
dna_adj<-unique(dna_adj)
dna_adj_mat<-as.matrix(dna_adj)
mode(dna_adj_mat) <- "numeric"

## Colors to distinguish between replicate pairs
subject_colors <- c(
  "pair1" = "#E41A1C", "pair2" = "#377EB8", "pair3" = "#4DAF4A", 
  "pair4" = "#984EA3", "pair5" = "#FF7F00", "pair6" = "#FFFF33", 
  "pair7" = "#A65628", "pair8" = "#F781BF", "pair9" = "#999999", 
  "pair10" = "#66C2A5", "pair11" = "#FC8D62", "pair12" = "#8DA0CB"
)

# Set up annotations for heat map
#View(colData)
annotation_col <- data.frame(condition = colData$condition, pair=colData$subject)
rownames(annotation_col) <- rownames(colData) 
head(annotation_col) #match replicate label to age

### LOG FOLD CHANGE ANNOTATION

dna_info<- selected_genes_adj[selected_genes_adj$Gene_Name %in% dna_repair, ]
dna_info$logBM<-log10(dna_info$baseMean)
dna_info<-unique(dna_info[,c("log2FoldChange", "logBM", "padj", "Gene_ID", "Gene_Name")])
rownames(dna_info) <- dna_info$Gene_Name

# Create a color palette for the LFC values
lfc_colors <- colorRampPalette(c("purple", "white", "darkgreen"))(100)
# bM_colors<- colorRampPalette(c("white", "black"))(100)


# Create a data frame for annotation, including LFC and BM
annotation_row <- data.frame(
  LFC = dna_info$log2FoldChange 
  #, LBM = logfc_info$logBM
)

# extract just the logFC info
rownames(annotation_row) <- dna_info$Gene_Name # Ensure row names match the heatmap data

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
all(rownames(annotation_row) == rownames(dna_adj_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(dna_adj_mat))


#-------------------------------------------------------------------------------
## Plot heatmap
### p < 0.1

#pdf("temp_figs/heatmap_DESEQadj_p<0.1.pdf", width = 8, height = 12)
#jpeg("figures/heatmap_DESEQadj_p<0.1.jpeg", width = 12, height = 18, units = "in", res = 300, quality = 85)

dna_map<- pheatmap(
 dna_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  legend= FALSE,
  annotation_row = annotation_row,
  #annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_legend = FALSE, 
  fontsize_row = 10, # Adjust the font size of row names
  fontsize_col = 10, # Adjust the font size of column names
  cellheight = 12)

## legend needs to be the whole thing...
legend_map<- pheatmap(
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

#dev.off()

#---------------------------------------------------------------------------
# ion regulation

# scale and convert to matrix
ion_adj<-t(scale(t(ion_reg))) 
#scale so each feature has the same mean/variance for visualization purposes 
# in the heatmap (prevents features with higher overall expression from washing 
# out signals of lower expression features)
ion_adj<-as.data.frame(ion_adj)
ion_adj<-unique(ion_adj)
ion_adj_mat<-as.matrix(ion_adj)
mode(ion_adj_mat) <- "numeric"

## Colors to distinguish between replicate pairs
# subject_colors <- c(
#   "pair1" = "#E41A1C", "pair2" = "#377EB8", "pair3" = "#4DAF4A", 
#   "pair4" = "#984EA3", "pair5" = "#FF7F00", "pair6" = "#FFFF33", 
#   "pair7" = "#A65628", "pair8" = "#F781BF", "pair9" = "#999999", 
#   "pair10" = "#66C2A5", "pair11" = "#FC8D62", "pair12" = "#8DA0CB"
# )

# Set up annotations for heat map
#View(colData)
# annotation_col <- data.frame(condition = colData$condition, pair=colData$subject)
# rownames(annotation_col) <- rownames(colData) 
# head(annotation_col) #match replicate label to age

### LOG FOLD CHANGE ANNOTATION

ion_info<- selected_genes_adj[selected_genes_adj$Gene_Name %in% ion, ]
# ion_info$logBM<-log10(ion_info$baseMean)
ion_info<-unique(ion_info[,c("log2FoldChange", "padj", "Gene_ID", "Gene_Name")])
rownames(ion_info) <- ion_info$Gene_Name


# Create a color palette for the LFC values
lfc_colors <- colorRampPalette(c("purple", "white", "darkgreen"))(100)
# bM_colors<- colorRampPalette(c("white", "black"))(100)


# Create a data frame for annotation, including LFC and BM
annotation_row <- data.frame(
  LFC = ion_info$log2FoldChange 
  #, LBM = logfc_info$logBM
)

# extract just the logFC info
rownames(annotation_row) <- ion_info$Gene_Name # Ensure row names match the heatmap data

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
all(rownames(annotation_row) == rownames(ion_adj_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(ion_adj_mat))


#-------------------------------------------------------------------------------
## Plot heatmap
### p < 0.1

#pdf("temp_figs/heatmap_DESEQadj_p<0.1.pdf", width = 8, height = 12)
#jpeg("figures/heatmap_DESEQadj_p<0.1.jpeg", width = 12, height = 18, units = "in", res = 300, quality = 85)

dna_map<-pheatmap(
  ion_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  legend = FALSE, 
  show_colnames = FALSE,
  annotation_row = annotation_row,
 # annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_legend = FALSE, 
  fontsize_row = 10, # Adjust the font size of row names
  fontsize_col = 10, # Adjust the font size of column names
  cellheight = 12)

#dev.off()

#---------------------------------------------------------------------------
# transcription and translational regulation

# scale and convert to matrix
trans_adj<-t(scale(t(trans_reg2)))

#scale so each feature has the same mean/variance for visualization purposes 
# in the heatmap (prevents features with higher overall expression from washing 
# out signals of lower expression features)
trans_adj<-as.data.frame(trans_adj)
trans_adj<-unique(trans_adj)
trans_adj_mat<-as.matrix(trans_adj)
mode(trans_adj_mat) <- "numeric"


## Colors to distinguish between replicate pairs
# subject_colors <- c(
#   "pair1" = "#E41A1C", "pair2" = "#377EB8", "pair3" = "#4DAF4A", 
#   "pair4" = "#984EA3", "pair5" = "#FF7F00", "pair6" = "#FFFF33", 
#   "pair7" = "#A65628", "pair8" = "#F781BF", "pair9" = "#999999", 
#   "pair10" = "#66C2A5", "pair11" = "#FC8D62", "pair12" = "#8DA0CB"
# )

# Set up annotations for heat map
#View(colData)
# annotation_col <- data.frame(condition = colData$condition, pair=colData$subject)
# rownames(annotation_col) <- rownames(colData) 
# head(annotation_col) #match replicate label to age

### LOG FOLD CHANGE ANNOTATION

trans_info<- selected_genes_adj[selected_genes_adj$Gene_Name %in% trans_reg, ]
# ion_info$logBM<-log10(ion_info$baseMean)
trans_info<-unique(trans_info[,c("log2FoldChange", "padj", "Gene_ID", "Gene_Name")])
rownames(trans_info) <- trans_info$Gene_Name


# Create a data frame for annotation, including LFC and BM
annotation_row <- data.frame(
  LFC = trans_info$log2FoldChange 
  #, LBM = logfc_info$logBM
)
# extract just the logFC info
rownames(annotation_row) <- trans_info$Gene_Name # Ensure row names match the heatmap data

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
all(rownames(annotation_row) == rownames(trans_adj_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(trans_adj_mat))

library(pheatmap)

trans_map<-pheatmap(
  trans_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  # annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_legend = FALSE, 
  legend = FALSE, 
  fontsize_row = 10, # Adjust the font size of row names
  fontsize_col = 10, # Adjust the font size of column names
  cellheight = 12)

#---------------------------------------------------------------------------
# biosynth

# scale and convert to matrix
bio_adj<-t(scale(t(biosynth2)))

#scale so each feature has the same mean/variance for visualization purposes 
# in the heatmap (prevents features with higher overall expression from washing 
# out signals of lower expression features)
bio_adj<-as.data.frame(bio_adj)
bio_adj<-unique(bio_adj)
bio_adj_mat<-as.matrix(bio_adj)
mode(bio_adj_mat) <- "numeric"


## Colors to distinguish between replicate pairs
# subject_colors <- c(
#   "pair1" = "#E41A1C", "pair2" = "#377EB8", "pair3" = "#4DAF4A", 
#   "pair4" = "#984EA3", "pair5" = "#FF7F00", "pair6" = "#FFFF33", 
#   "pair7" = "#A65628", "pair8" = "#F781BF", "pair9" = "#999999", 
#   "pair10" = "#66C2A5", "pair11" = "#FC8D62", "pair12" = "#8DA0CB"
# )

# Set up annotations for heat map
#View(colData)
# annotation_col <- data.frame(condition = colData$condition, pair=colData$subject)
# rownames(annotation_col) <- rownames(colData) 
# head(annotation_col) #match replicate label to age

### LOG FOLD CHANGE ANNOTATION
bio_info<- selected_genes_adj[selected_genes_adj$Gene_Name %in% biosynth, ]
# ion_info$logBM<-log10(ion_info$baseMean)
bio_info<-unique(bio_info[,c("log2FoldChange", "padj", "Gene_ID", "Gene_Name")])
rownames(bio_info) <- bio_info$Gene_Name

# Create a data frame for annotation, including LFC and BM
annotation_row <- data.frame(
  LFC = bio_info$log2FoldChange 
  #, LBM = logfc_info$logBM
)
# extract just the logFC info
rownames(annotation_row) <- bio_info$Gene_Name # Ensure row names match the heatmap data

### check that everything looks right
all(rownames(annotation_row) == rownames(bio_adj_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(bio_adj_mat))

library(pheatmap)

bio_map<-pheatmap(
  bio_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  # annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  legend = FALSE, 
  annotation_legend = FALSE, 
  fontsize_row = 10, # Adjust the font size of row names
  fontsize_col = 10, # Adjust the font size of column names
  cellheight = 12)

#---------------------------------------------------------------------------
# protein regulation

# scale and convert to matrix
prot_adj<-t(scale(t(protein_reg)))

#scale so each feature has the same mean/variance for visualization purposes 
# in the heatmap (prevents features with higher overall expression from washing 
# out signals of lower expression features)
prot_adj<-as.data.frame(prot_adj)
prot_adj<-unique(prot_adj)
prot_adj_mat<-as.matrix(prot_adj)
mode(prot_adj_mat) <- "numeric"


## Colors to distinguish between replicate pairs
# subject_colors <- c(
#   "pair1" = "#E41A1C", "pair2" = "#377EB8", "pair3" = "#4DAF4A", 
#   "pair4" = "#984EA3", "pair5" = "#FF7F00", "pair6" = "#FFFF33", 
#   "pair7" = "#A65628", "pair8" = "#F781BF", "pair9" = "#999999", 
#   "pair10" = "#66C2A5", "pair11" = "#FC8D62", "pair12" = "#8DA0CB"
# )

# Set up annotations for heat map
#View(colData)
# annotation_col <- data.frame(condition = colData$condition, pair=colData$subject)
# rownames(annotation_col) <- rownames(colData) 
# head(annotation_col) #match replicate label to age

### LOG FOLD CHANGE ANNOTATION
prot_info<- selected_genes_adj[selected_genes_adj$Gene_Name %in% protein, ]
# ion_info$logBM<-log10(ion_info$baseMean)
prot_info<-unique(prot_info[,c("log2FoldChange", "padj", "Gene_ID", "Gene_Name")])
rownames(prot_info) <- prot_info$Gene_Name

# Create a data frame for annotation, including LFC and BM
annotation_row <- data.frame(
  LFC = prot_info$log2FoldChange 
  #, LBM = logfc_info$logBM
)
# extract just the logFC info
rownames(annotation_row) <- prot_info$Gene_Name # Ensure row names match the heatmap data

### check that everything looks right
all(rownames(annotation_row) == rownames(prot_adj_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(prot_adj_mat))

library(pheatmap)

prot_map<-pheatmap(
  prot_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  # annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_legend = FALSE, 
  fontsize_row = 10, # Adjust the font size of row names
  fontsize_col = 10, # Adjust the font size of column names
  legend = FALSE, 
  cellheight = 12)

#---------------------------------------------------------------------------
# structural 

# scale and convert to matrix
str_adj<-t(scale(t(structural)))

#scale so each feature has the same mean/variance for visualization purposes 
# in the heatmap (prevents features with higher overall expression from washing 
# out signals of lower expression features)
str_adj<-as.data.frame(str_adj)
str_adj<-unique(str_adj)
str_adj_mat<-as.matrix(str_adj)
mode(str_adj_mat) <- "numeric"


## Colors to distinguish between replicate pairs
# subject_colors <- c(
#   "pair1" = "#E41A1C", "pair2" = "#377EB8", "pair3" = "#4DAF4A", 
#   "pair4" = "#984EA3", "pair5" = "#FF7F00", "pair6" = "#FFFF33", 
#   "pair7" = "#A65628", "pair8" = "#F781BF", "pair9" = "#999999", 
#   "pair10" = "#66C2A5", "pair11" = "#FC8D62", "pair12" = "#8DA0CB"
# )

# Set up annotations for heat map
#View(colData)
# annotation_col <- data.frame(condition = colData$condition, pair=colData$subject)
# rownames(annotation_col) <- rownames(colData) 
# head(annotation_col) #match replicate label to age

### LOG FOLD CHANGE ANNOTATION
str_info<- selected_genes_adj[selected_genes_adj$Gene_Name %in% structure, ]
# ion_info$logBM<-log10(ion_info$baseMean)
str_info<-unique(str_info[,c("log2FoldChange", "padj", "Gene_ID", "Gene_Name")])
rownames(str_info) <- str_info$Gene_Name

# Create a data frame for annotation, including LFC and BM
annotation_row <- data.frame(
  LFC = str_info$log2FoldChange 
  #, LBM = logfc_info$logBM
)
# extract just the logFC info
rownames(annotation_row) <- str_info$Gene_Name # Ensure row names match the heatmap data

### check that everything looks right
all(rownames(annotation_row) == rownames(str_adj_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(str_adj_mat))


str_map<-pheatmap(
  str_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  # annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  legend = FALSE, 
  annotation_legend = FALSE, 
  fontsize_row = 10, # Adjust the font size of row names
  fontsize_col = 10, # Adjust the font size of column names
  cellheight = 12)

#---------------------------------------------------------------------------
# other

# scale and convert to matrix
other_adj<-t(scale(t(other2)))

#scale so each feature has the same mean/variance for visualization purposes 
# in the heatmap (prevents features with higher overall expression from washing 
# out signals of lower expression features)
other_adj<-as.data.frame(other_adj)
other_adj<-unique(other_adj)
other_adj_mat<-as.matrix(other_adj)
mode(other_adj_mat) <- "numeric"
View(other_adj_mat)

## Colors to distinguish between replicate pairs
# subject_colors <- c(
#   "pair1" = "#E41A1C", "pair2" = "#377EB8", "pair3" = "#4DAF4A", 
#   "pair4" = "#984EA3", "pair5" = "#FF7F00", "pair6" = "#FFFF33", 
#   "pair7" = "#A65628", "pair8" = "#F781BF", "pair9" = "#999999", 
#   "pair10" = "#66C2A5", "pair11" = "#FC8D62", "pair12" = "#8DA0CB"
# )

# Set up annotations for heat map
#View(colData)
# annotation_col <- data.frame(condition = colData$condition, pair=colData$subject)
# rownames(annotation_col) <- rownames(colData) 
# head(annotation_col) #match replicate label to age

### LOG FOLD CHANGE ANNOTATION
other_info<- selected_genes_adj[selected_genes_adj$Gene_Name %in% other, ]
# ion_info$logBM<-log10(ion_info$baseMean)
other_info<-unique(other_info[,c("log2FoldChange", "padj", "Gene_ID", "Gene_Name")])
rownames(other_info) <- other_info$Gene_Name

# Create a data frame for annotation, including LFC and BM
annotation_row <- data.frame(
  LFC = other_info$log2FoldChange 
  #, LBM = logfc_info$logBM
)
# extract just the logFC info
rownames(annotation_row) <- other_info$Gene_Name # Ensure row names match the heatmap data

### check that everything looks right
all(rownames(annotation_row) == rownames(other_adj_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(other_adj_mat))


other_map<-pheatmap(
  other_adj_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  # annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_legend = FALSE, 
  fontsize_row = 10, # Adjust the font size of row names
  fontsize_col = 10, # Adjust the font size of column names
  legend = FALSE, 
  cellheight = 12)

#-------------------------------------------------------------------------------
### Convert pheatmaps to complex heatmaps


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)
library(dplyr)

ComplexHeatmap::pheatmap(other_map)



