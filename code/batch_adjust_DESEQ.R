### Correct for batch effect and create DESEQ object

### correct for batch effect from 2 different sorting groups

## convert back to dataframe after batch correction to level data correctly

#Read in Gene Count Matrix, subset
counts <- read.table("data/gene_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
CountData<-counts[,-(1:5)]
CountData<-as.matrix(CountData) ## this is important--combat seq doesn't like dataframes

## read in design file and designate as factors/matrix
colData <- read.table("data/design.txt", header = TRUE, sep = "\t", row.names = 1)
colData<-as.matrix(colData)

# batch information
batch <- c(rep(1, 6), rep(2, 6), rep(1, 6), rep(2, 6)) #Organize data into batch 1 (pairs 1-6, sorted on Day 1) and batch 2 (pairs 7-12, sorted on Day 2)
age<-c(rep(1, 12), rep(2, 12)) # indicate whether the replicate is "old" (1) or "young" (2)
#pair<-rep(1:12, 2) #indicate the paired nature of the old and young replicates # don't include this in CombatSeq

# Run ComBat_seq
adj_counts <- ComBat_seq(CountData, batch=batch, group=age) #ignore the covariates...package is not super clear on what they are used for, and we don't want to include pair in both ComBat-seq and the DEseq model

write.table(adj_counts, file="temp/gcm_combatseq.txt")

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
any(is.na(dds_adj)) # check for NAs, there are none so we can skip the next line
#res_adj<-results(dds_adj, cooksCutoff = FALSE) #save results table #prevents cooks cuttoff from assigning NA values (can also test independentFiltering to false if still having issues).
res_df<-as.data.frame(res_adj) 
View(res_df)

### mean basemean of p-values below 0.05
# Filter for rows with p-value > 0.05
filtered_res <- res_df[res_df$padj < 0.05, ]

# Calculate the mean of the baseMean values
mean_baseMean <- mean(filtered_res$baseMean, na.rm = TRUE)

# Print the result
mean_baseMean

head(results(dds_adj, tidy=TRUE)) #cleaner results table
summary(res_adj) #summary table
res_adj <- res_adj[order(res_adj$padj),] #sort by p-value

#reorder
resOrdered_adj <- res_adj[order(res_adj$pvalue),]
head(resOrdered_adj)

write.csv(as.data.frame(resOrdered_adj), 
          file="temp/rnaseq_results_batch_adjusted.csv") ### this contains 
## results from the dds object

### additional options for plotting and visualization in Analysis_eNotebook_Part2_DGE.rmd file
