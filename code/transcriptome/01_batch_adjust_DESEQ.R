### Correct for batch effect and create DESEQ object
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## libraries
#-------------------------------------------------------------------------------

#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("sva")
library(sva)
library(DESeq2)
#-------------------------------------------------------------------------------

### correct for batch effect from 2 different sorting groups
## convert back to dataframe after batch correction to level data correctly

#-------------------------------------------------------------------------------
#Read in Gene Count Matrix, subset

counts <- read.table("data/raw/gene_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
CountData<-counts[,-(1:5)]
CountData<-as.matrix(CountData) ## this is important--combat seq doesn't like dataframes

## read in design file and designate as factors/matrix
colData <- read.table("data/design_files/design.txt", header = TRUE, sep = "\t", row.names = 1)
colData$batch <- c(
  "B1", "B1", "B1", "B1","B1", "B1", 
  "B2", "B2", "B2", "B2", "B2", "B2",
  "B1", "B1", "B1", "B1","B1", "B1", 
  "B2", "B2", "B2", "B2", "B2", "B2"
)
colData<-as.matrix(colData)


#-------------------------------------------------------------------------------
# batch information
batch <- c(rep(1, 6), rep(2, 6), rep(1, 6), rep(2, 6)) #Organize data into batch 1 (pairs 1-6, sorted on Day 1) and batch 2 (pairs 7-12, sorted on Day 2)
age<-c(rep(1, 12), rep(2, 12)) # indicate whether the replicate is "old" (1) or "young" (2)
#pair<-rep(1:12, 2) #indicate the paired nature of the old and young replicates # don't include this in CombatSeq
#-------------------------------------------------------------------------------

# Run ComBat_seq
adj_counts <- ComBat_seq(CountData, batch=batch, group=age) #ignore the covariates...package is not super clear on what they are used for, and we don't want to include pair in both ComBat-seq and the DEseq model
View(adj_counts)
#-------------------------------------------------------------------------------
#write
write.table(adj_counts, file="temp/transcriptome/gcm_combatseq.txt")
#-------------------------------------------------------------------------------

### Create DESEQ object
colData2<-as.data.frame(colData)
colData2$condition <- factor(colData2$condition)
levels(colData2$condition) #
colData2$condition <- relevel(colData2$condition, ref = "young") # we want to treat the young as the reference and discuss changes in old
dds_adj<-DESeqDataSetFromMatrix(countData= adj_counts, colData=colData2, design= ~subject + condition)


#remove genes that have low mapping (>=5), >=13 makes sure that genes that have near 0 expression in only one treatment aren't excluded, since these are genes of interest
keep<- rowSums(counts(dds_adj) >=5)>=13 #and then number of samples that have >=13 #go check the manual 
dds_adj<- dds_adj[keep,]

nrow(dds_adj) #filtering removed about 1000 genes #5620 genes left

### Run DESeq analysis
dds_adj<-DESeq(dds_adj) 

### need to normalize after DESEQ for visualization (not analysis)
normalized_counts_adj <- counts(dds_adj, normalized = TRUE) # Get normalized counts

any(is.na(dds_adj)) # check for NAs, there are none so we can skip the next line
res_adj<-results(dds_adj, cooksCutoff = FALSE) #save results table #prevents cooks cuttoff from assigning NA values (can also test independentFiltering to false if still having issues).
write.table(res_adj, file="temp/transcriptome/dds_object_adjusted.txt")
dds_adj1<-as.data.frame(res_adj)
### mean basemean of p-values below 0.05
# Filter for rows with p-value > 0.05
filtered_res <- dds_adj1[dds_adj1$padj < 0.05, ]
View(filtered_res)

# Calculate the mean of the baseMean values
mean_baseMean <- mean(filtered_res$baseMean, na.rm = TRUE)

# Print the result
mean_baseMean

summary(dds_adj1) #summary table
dds_adj1 <- dds_adj1[order(dds_adj1$padj),] #sort by p-value

#reorder
resOrdered_adj <- dds_adj1[order(dds_adj1$pvalue),]
head(resOrdered_adj)
nrow(resOrdered_adj)

#-------------------------------------------------------------------------------
# write data
write.csv(as.data.frame(resOrdered_adj), 
          file="data/clean/rnaseq_results_batch_adjusted.csv") ### this contains 
## results from the dds object

write.csv(normalized_counts_adj, file="data/clean/normalized_counts_deseq.csv")

### additional options for plotting and visualization in Analysis_eNotebook_Part2_DGE.rmd file


#--------------- PCAs

# VISUALIZATION ONLY-- transform  data to make it homoskedastic (variance of the residual is constant)  #this is JUST for visualization
rlog_dds<-rlog(dds_adj)
df <-as.data.frame(assay(rlog_dds)[,12:13])

test<-ggplot(df, aes(x= RNAseq_OLD_rep12.bam, y=RNAseq_YOUNG_rep01.bam)) +geom_hex(bins=100, colour="orange", fill="black")+coord_fixed()+theme_classic()
ggsave("figures/test_adjusted.pdf", test)

#PCA
pca_plot<-plotPCA(rlog_dds)
pca2<-plotPCA(rlog_dds, intgroup = "batch")
ggsave("figures/pca_RNAseq_adjusted.pdf", pca_plot)


pcaData <- plotPCA(rlog_dds,
                   intgroup = c("batch","condition"),
                   returnData = TRUE)

p2<-ggplot(pcaData,
          aes(PC1, PC2,
              color = batch, 
              shape = condition)) +
  geom_point(size = 4)

p2<-p2+theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 14)) 



ggsave("figures/pca_RNAseq_batch_adjusted.pdf", p2)

?plotPCA

#generate a loading plot
#transpose data
t_rld<-t(assay((rlog_dds)))
t_rld[1:3,1:4]

pca<-prcomp(t_rld)
plot(pca)

#loadings
names(pca)
loadings<-as.data.frame(pca$rotation)
head(loadings)
loadingplot(loadings$PC1, threshold = 0.1)
plot_load<-loadingplot(loadings$PC1,threshold= 0.1)
row.names(dds)[plot_load$var.names] ###BUT...the PCA doesn't look awesome, so I'm not sure how useful this really is



