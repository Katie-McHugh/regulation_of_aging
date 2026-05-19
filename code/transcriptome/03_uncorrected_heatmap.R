### Correct for batch effect and create DESEQ object
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## libraries
#-------------------------------------------------------------------------------

#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("sva")
library(sva)
library(DESeq2)
#-------------------------------------------------------------------------------

### did not correct for batch effect from 2 different sorting groups
## convert back to dataframe after batch correction to level data correctly

#-------------------------------------------------------------------------------
#Read in Gene Count Matrix, subset

counts <- read.table("data/raw/gene_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
CountData<-counts[,-(1:5)]
CountData<-as.matrix(CountData) ## this is important--combat seq doesn't like dataframes
View(CountData)

## read in design file and designate as factors/matrix
colData <- read.table("data/design_files/design.txt", header = TRUE, sep = "\t", row.names = 1)
colData$batch <- c(
  "Batch 1", "Batch 1", "Batch 1", "Batch 1", "Batch 1", "Batch 1", 
  "Batch 2", "Batch 2", "Batch 2", "Batch 2", "Batch 2", "Batch 2",
  "Batch 1", "Batch 1", "Batch 1", "Batch 1", "Batch 1", "Batch 1", 
  "Batch 2", "Batch 2", "Batch 2", "Batch 2", "Batch 2", "Batch 2"
)
colData<-as.matrix(colData)
View(colData)

### Create DESEQ object
colData2<-as.data.frame(colData)
colData2$condition <- factor(colData2$condition)
levels(colData2$condition) 
colData2$condition <- relevel(colData2$condition, ref = "young") # we want to treat the young as the reference and discuss changes in old
dds<-DESeqDataSetFromMatrix(countData= CountData, colData=colData2, design= ~subject + condition)


#remove genes that have low mapping (>=5), >=13 makes sure that genes that have near 0 expression in only one treatment aren't excluded, since these are genes of interest
keep<- rowSums(counts(dds) >=5)>=13 #and then number of samples that have >=13 #go check the manual 
dds<- dds[keep,]

nrow(dds) #filtering removed about 1000 genes #5620 genes left

### Run DESeq analysis
dds<-DESeq(dds) 

### need to normalize after DESEQ for visualization (not analysis)
normalized_counts <- counts(dds, normalized = TRUE) # Get normalized counts
any(is.na(dds)) # check for NAs, there are none so we can skip the next line
res<-results(dds, cooksCutoff = FALSE) #save results table #prevents cooks cuttoff from assigning NA values (can also test independentFiltering to false if still having issues).

dds_1<-as.data.frame(res)

### mean basemean of p-values below 0.05
# Filter for rows with p-value > 0.05
filtered_res <- dds_1[dds_1$padj < 0.05, ]
filtered_res2 <- dds_1[dds_1$padj < 0.1, ]

# Calculate the mean of the baseMean values
mean_baseMean <- mean(filtered_res$baseMean, na.rm = TRUE)

# Print the result
mean_baseMean

summary(dds_1) #summary table
dds_1 <- dds_1[order(dds_1$padj),] #sort by p-value

#reorder
resOrdered <- dds_1[order(dds_1$pvalue),]
head(resOrdered)
nrow(resOrdered)
sigs_res<- subset(resOrdered, padj<0.1)
nrow(sigs_res)

#-------------------------------------------------------------------------------
# write data
#write.csv(as.data.frame(resOrdered), 
         # file="data/clean/rnaseq_results_batch_adjusted_uncorrected.csv") ### this contains 
## results from the dds object

#write.csv(normalized_counts, file="data/clean/normalized_counts_deseq_uncorrected.csv")

### additional options for plotting and visualization in Analysis_eNotebook_Part2_DGE.rmd file

#-------------------------------------------------------------------------------
## heatmap for uncorrected data

library(pheatmap)
library(tidyverse)
# 
# okabe_ito_palette <- c("#D55E00",
#                        "#F0E442",
#                        "#56B4E9",
#                        "#0072B2", 
#                        "#E69F00")
# "#009E73",
# "#CC79A7"

## Load and Organize data

### Load normalized counts for visualization
norm_dds<-as.data.frame(normalized_counts)
norm_dds <- rownames_to_column(norm_dds, var = "Gene_ID")



### Load list of genes for heatmap
selected_genes<-sigs_res
nrow(selected_genes)
selected_genes <- rownames_to_column(selected_genes, var = "Gene_ID")

# change "old" to "aged"
colData<-as.data.frame(colData)
colData <- colData %>%
  mutate(across(where(is.character), ~ gsub("old", "aged", .)))
colData$condition <- factor(colData$condition)
colData$condition <- relevel(colData$condition, ref = "young") # we want to treat the young as the reference and discuss changes in old


### Load in gene key
key<-read.table("temp/transcriptome/key_geneIDtoName.txt")

### subset norm_dds to significant genes
sel<-as.data.frame(selected_genes$Gene_ID)
colnames(sel)[colnames(sel) == "selected_genes$Gene_ID"] <- "Gene_ID"
norm_sigs<-merge(sel, norm_dds, by="Gene_ID", all.x=TRUE)
View(norm_sigs)

# Set key for old vs young replicates
age_colors<- c( "aged" = "black", "young"= "grey")
age_cols<-list(condition = age_colors)

### Convert gene ID to gene name, and make that the row name
norm_sigs2<-merge(norm_sigs, key, by="Gene_ID", all.x = TRUE)
View(norm_sigs2)
## rename using SGD
norm_sigs2$Gene_Name[norm_sigs2$Gene_ID == "YLR050C"] <- "EMA19" ## Name recognized in SGD
norm_sigs2$Gene_Name[norm_sigs2$Gene_ID == "YCR015C"] <- "CTO1" ## Name recognized in SGD

rownames(norm_sigs2) <- norm_sigs2$Gene_Name

# Remove the non-numeric columns
norm_sigs2$Gene_ID <- NULL
norm_sigs2$Gene_Name <- NULL

# scale and convert to matrix
norm<-t(scale(t(norm_sigs2))) 
#scale so each feature has the same mean/variance for visualization purposes 
# in the heatmap (prevents features with higher overall expression from washing 
# out signals of lower expression features)
norm<-as.data.frame(norm)
norm<-unique(norm)
norm_mat<-as.matrix(norm)
mode(norm_mat) <- "numeric"



## Colors to distinguish between replicate pairs
subject_colors <- c(
  "pair1" = "#E41A1C", "pair2" = "#377EB8", "pair3" = "#4DAF4A", 
  "pair4" = "#984EA3", "pair5" = "#FF7F00", "pair6" = "#FFFF33", 
  "pair7" = "#A65628", "pair8" = "#F781BF", "pair9" = "#999999", 
  "pair10" = "#66C2A5", "pair11" = "#FC8D62", "pair12" = "#8DA0CB"
)

# Set up annotations for heat map
annotation_col <- data.frame(condition = colData$condition, pair=colData$subject)
rownames(annotation_col) <- rownames(colData) 
head(annotation_col) #match replicate label to age


### LOG FOLD CHANGE ANNOTATION
selected_genes$logBM<-log10(selected_genes$baseMean)

selected_genes2<-merge(selected_genes, key, by="Gene_ID", all.x = TRUE)
head(selected_genes2)
selected_genes2$Gene_Name[selected_genes2$Gene_ID == "YLR050C"] <- "EMA19" ## Name recognized in SGD
selected_genes2$Gene_Name[selected_genes2$Gene_ID == "YCR015C"] <- "CTO1" ## Name recognized in SGD

logfc_info<-unique(selected_genes2[,c("log2FoldChange", "logBM", "padj", "Gene_ID", "Gene_Name")])
rownames(logfc_info) <- logfc_info$Gene_Name



# Create a color palette for the LFC values
lfc_colors <- colorRampPalette(c("#CC79A7", "white", "#009E73"))(100)
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
all(rownames(annotation_row) == rownames(norm_mat))  # Check if row names match
all(rownames(annotation_col) == colnames(norm_mat))

View(norm_mat)
nrow(norm_mat)
nrow(annotation_row)


View(logfc_info)

#-------------------------------------------------------------------------------
## Plot heatmap
### p < 0.1

#pdf("figures/SuppFig3_heatmap.pdf", width = 8, height = 12)
tiff(file = "figures/Supp_FIG_heatmap_UNCORRECTED.tiff", height = 12, width = 12,
     units = "in",
     res = 400,
     compression = "lzw")

pheatmap(
  norm_mat, # Scale the data by rows (genes)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  treeheight_row = 20,  
  annotation_legend = TRUE, 
  fontsize = 10,           # Increase overall font size
  fontsize_row = 10,       # Row name font size
  fontsize_col = 10,       # Column name font size
  cellheight = 10,         # Increase cell height
  cellwidth = 20         # Increase cell width if needed
)

dev.off()

#--------------- PCAs for dds object -------------------#

# VISUALIZATION ONLY-- transform  data to make it homoskedastic (variance of the residual is constant)  #this is JUST for visualization
rlog_dds<-rlog(dds)

## test transformation #better than vst
df <-as.data.frame(assay(rlog_dds)[,12:13])
test<-ggplot(df, aes(x= RNAseq_OLD_rep12.bam, y=RNAseq_YOUNG_rep01.bam)) +geom_hex(bins=100, colour="orange", fill="black")+coord_fixed()+theme_classic()
ggsave("figures/test_adjusted.pdf", test)

#PCA

#pca_dds<-plotPCA(rlog_dds, intgroup = "batch")
#pca_dds2<-plotPCA(rlog_dds, intgroup = "condition")

## ggplot is easier to work with

pcaData_dds <- plotPCA(rlog_dds,
                       intgroup = c("batch","condition"),
                       returnData = TRUE)

#pcaData_dds$batch <- ifelse(pcaData_dds$batch == "B1", "Batch 1", "Batch 2")

head(pcaData_dds, n=12)
attr(pcaData_dds, "percentVar")

percentVar <- round(100 * attr(pcaData_dds, "percentVar"), 2)

View(pcaData_dds)
pca_dds1 <- ggplot(pcaData_dds, aes(PC1, PC2, color = batch)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Batch 1" = "steelblue", "Batch 2" = "tomato")) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(x = paste0('PC1: ', percentVar[1], '%'),
       y = paste0('PC2: ', percentVar[2], '%'),
       color = "Batch") +
  theme_bw()

pca_dds1


pca_dds2 <- ggplot(pcaData_dds, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("old" = "orange", "young" = "purple"),
                     breaks = c("old", "young"))+
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(x = paste0('PC1: ', percentVar[1], '%'),
       y = paste0('PC2: ', percentVar[2], '%'),
       color = "condition") +
  theme_bw()

pca_dds2

library(patchwork)

pca_dds_combined<-(pca_dds1 / pca_dds2 + plot_annotation(tag_levels = list(c('A', 'B'))))
pca_dds_combined
ggsave("temp/transcriptome/pca_dds_uncorrected.pdf", pca_dds_combined, width = 4.085, height = 8.06)

