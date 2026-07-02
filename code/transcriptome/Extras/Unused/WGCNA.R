#### WGCNA
## helpful tutorial: https://www.youtube.com/watch?v=
### for paired data: https://www.nature.com/articles/s41598-017-18705-z


install.packages('BiocManager')
library(BiocManager)
BiocManager::install('WGCNA')
BiocManager::install('flashClust')
BiocManager::install("GEOquery")
BiocManager::install("CorLevelPlot")
install.packages('GEOquery')

library(WGCNA)
library(flashClust)
library(curl)
library(GEOquery)
library(gridExtra)
library(tidyverse)
library(CorLevelPlot)

# load in data
data<-read.table("data/raw/gene_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
CountData<-data[,-(1:5)]
colnames(CountData) <- tolower(gsub("RNAseq_(.+)_rep([0-9]+)\\.bam", "\\1 \\2", colnames(CountData)))

gsg<-goodSamplesGenes(t(CountData))
summary(gsg)
gsg$allOK #FALSE = outliers are present

table(gsg$goodGenes) #266 outliers
table(gsg$goodSamples)

data <- CountData[gsg$goodGenes ==TRUE,]

htree<-hclust(dist(t(CountData)), method = "average")
plot(htree) ## shows outlier samples

pca<-prcomp(t(data))
pca.data<-pca$x
pca.data<-as.data.frame(pca.data)
pca.data$group <- gsub(" [0-9]+", "", rownames(pca.data))
pca.data$rep <- as.numeric(sub(".* ", "", rownames(pca.data)))
pca.data$batch <- ifelse(pca.data$rep <= 6, "Batch1", "Batch2")
View(pca.data)

pca.var<-pca$sdev^2
pca.var.percent<-round(pca.var/sum(pca.var)*100, digits=2)

pca1<-ggplot(pca.data, aes(PC1, PC2, color= batch)) +
  geom_point(size = 3) +
  # geom_text(label = rownames(pca.data), 
  #           vjust = -0.8,    # moves text up (negative) or down (positive)
  #           hjust = 0.5,     # horizontal centering
  #           size = 3) +      # smaller text (default is 3.88)) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +  # 15% padding on x axis
  scale_y_continuous(expand = expansion(mult = 0.15)) +  # 15% padding on y axis
  scale_color_manual(values = c("Batch1" = "steelblue", "Batch2" = "tomato")) +
  labs(x=paste0('PC1: ', pca.var.percent[1], '%'), 
       y=paste0('PC2: ', pca.var.percent[2], '%')) +
  theme_bw()

ggsave("temp/transcriptome/wgcna_pca_uncorrected.pdf", plot= pca1)
pca1


###########################
#batch adjusted
data2<-read.table(file="temp/transcriptome/gcm_combatseq.txt", row.names = 1) ### this contains 
#norm<-read.csv(file="data/clean/normalized_counts_deseq.csv", row.names=1)
colnames(data2) <- tolower(gsub("RNAseq_(.+)_rep([0-9]+)\\.bam", "\\1 \\2", colnames(data2)))

sum(data2== 0)
View(data2)

gsg2<-goodSamplesGenes(t(data2))
summary(gsg2)
gsg2$allOK #FALSE = outliers are present

table(gsg2$goodGenes) #266 outliers
table(gsg2$goodSamples)

data2 <- data2[gsg2$goodGenes ==TRUE,]
dev.off()
htree2<-hclust(dist(t(data2)), method = "average")
plot(htree2) ## shows outlier samples

pdf("temp/transcriptome/hierarchical_clustering_adj.pdf", width = 8, height = 5)

plot(htree2)

dev.off()

pca2<-prcomp(t(data2))
pca.data2<-pca2$x
pca.data2<-as.data.frame(pca.data2)
View(pca.data)
pca.data2$group <- gsub(" [0-9]+", "", rownames(pca.data2))
pca.data2$rep <- as.numeric(sub(".* ", "", rownames(pca.data2)))
pca.data2$batch <- ifelse(pca.data2$rep <= 6, "Batch1", "Batch2")

pca.var2<-pca2$sdev^2
pca.var.percent2<-round(pca.var2/sum(pca.var2)*100, digits=2)
View(pca.var.percent)

pca2<-ggplot(pca.data2, aes(PC1, PC2, color = batch)) +
  geom_point(size = 3) +
  #geom_text(label = rownames(pca.data2), 
  #          vjust = -0.8,    # moves text up (negative) or down (positive)
   #         hjust = 0.5,     # horizontal centering
    #        size = 3) +      # smaller text (default is 3.88)) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +  # 15% padding on x axis
  scale_y_continuous(expand = expansion(mult = 0.15)) +  # 15% padding on y axis
  scale_color_manual(values = c("Batch1" = "steelblue", "Batch2" = "tomato")) +
  labs(x=paste0('PC1: ', pca.var.percent2[1], '%'), 
       y=paste0('PC2: ', pca.var.percent2[2], '%'))+
  theme_bw()

pca2

ggsave("temp/transcriptome/wgcna_pca_batch_corrected.pdf", plot= pca2)
## 


### now color by treatment: 

pca3<-ggplot(pca.data, aes(PC1, PC2, color= group)) +
  geom_point(size = 3) +
  # geom_text(label = rownames(pca.data), 
  #           vjust = -0.8,    # moves text up (negative) or down (positive)
  #           hjust = 0.5,     # horizontal centering
  #           size = 3) +      # smaller text (default is 3.88)) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +  # 15% padding on x axis
  scale_y_continuous(expand = expansion(mult = 0.15)) +  # 15% padding on y axis
  scale_color_manual(values = c("old" = "orange", "young" = "purple")) +
  labs(x=paste0('PC1: ', pca.var.percent[1], '%'), 
       y=paste0('PC2: ', pca.var.percent[2], '%')) +
  theme_bw()

#ggsave("temp/transcriptome/wgcna_pca_uncorrected_agecolor.pdf", plot= pca1)
pca3

pca4<-ggplot(pca.data2, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  #geom_text(label = rownames(pca.data2), 
  #          vjust = -0.8,    # moves text up (negative) or down (positive)
  #         hjust = 0.5,     # horizontal centering
  #        size = 3) +      # smaller text (default is 3.88)) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +  # 15% padding on x axis
  scale_y_continuous(expand = expansion(mult = 0.15)) +  # 15% padding on y axis
  scale_color_manual(values = c("old" = "orange", "young" = "purple")) +
  labs(x=paste0('PC1: ', pca.var.percent2[1], '%'), 
       y=paste0('PC2: ', pca.var.percent2[2], '%'))+
  theme_bw()

pca4

library(patchwork)

pca_combined<-(pca1 / pca3 | pca2 / pca4) + plot_annotation(tag_levels = list(c('E', 'F', 'G', 'H')))
pca_combined
ggsave("temp/transcriptome/wgcna_pca.pdf", plot= pca_combined)

# get metadata
phenodata<-read.table("data/design_files/design.txt", header=TRUE, row.names = "sample")

head(phenodata)
head(data)

## read in design file and designate as factors/matrix
coldata2 <- read.table("data/design_files/design.txt", header = TRUE, sep = "\t", row.names = 1)
coldata1$batch <- c(
  "Batch 1", "Batch 1", "Batch 1", "Batch 1", "Batch 1", "Batch 1", 
  "Batch 2", "Batch 2", "Batch 2", "Batch 2", "Batch 2", "Batch 2",
  "Batch 1", "Batch 1", "Batch 1", "Batch 1", "Batch 1", "Batch 1", 
  "Batch 2", "Batch 2", "Batch 2", "Batch 2", "Batch 2", "Batch 2"
)




## remove outliers - old 2 and 7, along with paired young samples.
exclude<-c('old 01', 'old 07', 'young 01', 'young 07')
subset<-data[,!(colnames(data) %in% exclude)]

head(coldata2)

gsub("RNAseq_(.+)_rep[0-9]+\\.bam", "\\1", rownames(coldata2))
#rename
coldata2<-as.data.frame(coldata2)
rownames(coldata2) <- tolower(gsub("RNAseq_(.+)_rep([0-9]+)\\.bam", "\\1 \\2", rownames(coldata2)))

coldata2<-coldata2[!(rownames(coldata2) %in% exclude),]
head(coldata2)
head(subset)

all(rownames(coldata2) %in% colnames(subset))

## dds

dds<-DESeqDataSetFromMatrix(countData = subset,
                            colData = coldata2,
                            design = ~1)


dds_filtered <-dds[rowSums(counts(dds) >=5) >=13,]
dds_norm<-rlog(dds_filtered) #vs vst

norm.counts <- assay(dds_norm) %>%
  t()

## this went fast, need to come back and see what is happening
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
sft<-pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)

sft.data <-sft$fitIndices
View(sft.data) ## we want maximum r.sq and minimum connectivity

## visualize

a1<- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) + 
         geom_point() + 
         geom_text(nudge_y = 0.1) + 
         geom_hline(yintercept =0.8, color = 'red') + 
         labs(x = 'Power', y= 'Scale free topology model fit, signed R^2') +
         theme_bw()
a1

a2<- ggplot(sft.data, aes(Power, mean.k., label = Power)) + 
  geom_point() + 
  geom_text(nudge_y = 0.1) + 
  labs(x = 'Power', y= 'Mean connectivity') +
  theme_bw()

library(gridExtra)
grid.arrange(a1, a2, nrow = 2) ## pick one above the red line, but keep connectivity low

## now we've selected our power 

norm.counts[] <-sapply(norm.counts, as.numeric)
soft_power<- 16
temp_cor<-cor

## memory estimate
cor <-WGCNA::cor
bwnet<-blockwiseModules(norm.counts, 
                 maxBlockSize = 14000, 
                 TOMType = "signed", 
                 power = soft_power,
                 mergeCutHeght = 0.25, 
                 numericLabels = FALSE, 
                 randomSeed = 1234,
                 verbose = 3)

cor <-temp_cor

## s. module eigengenes
module_eigengenes<-bwnet$MEs
head(module_eigengenes)

##genes per module 
table(bwnet$colors)

## plot
plotDendroAndColors(bwnet$dendrograms, cbind(bwnet$unmergedColors, bwnet$colors), 
                    c("unmerged", "merged"),
                    dendroLabels = FALSE, 
                    addGuide = TRUE,
                    hang = 0.03, 
                    guideHang = 0.05)


## on to part 2 of the tutorial: 


