### Generating DGE Heatmap

### Load normalized counts for visualization
norm_dds<-read.csv("temp/transcriptome/normalized_counts_deseq.csv")

### Load list of genes for headmap
selected_genes_adj<-read.csv("temp/transcriptome/RNA_genes_p<0.1.csv", row.names= "X")
selected_genes_adj2<-read.csv("temp/transcriptome/RNA_genes_p<0.05.csv", row.names= "X")

### load in design file
colData<-read.table("data/design.txt", header=TRUE, row.names = "sample")

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
norm_adj<-t(scale(t(norm_sigs2))) #scale so each feature has the same mean/variance for visualization purposes in the heatmap (prevents features with higher overall expression from washing out signals of lower expression features)
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
  LFC = lfc_colors#,  # Use your color palette for LFC
  # LBM = bM_colors
)

### check that everything looks right
all(rownames(annotation_row) == rownames(norm_adj_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(norm_adj_mat))

### p < 0.1

#pdf("temp_figs/heatmap_DESEQadj_p<0.1.pdf", width = 8, height = 12)
jpeg("figures/heatmap_DESEQadj_p<0.01.jpeg", width = 12, height = 18, units = "in", res = 300, quality = 85)

pheatmap(
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

dev.off()
