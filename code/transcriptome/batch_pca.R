### PCAs of corrected and uncorrected data

## load in data
dds<-read.table(file="temp/transcriptome/dds_object_uncorrected.txt")
dds_adj<-read.table(file="temp/transcriptome/dds_object_adjusted.txt")

library(DESeq2)

# VISUALIZATION ONLY-- transform  data to make it homoskedastic (variance of the residual is constant)  #this is JUST for visualization
rlog_dds<-rlog(dds) #rlog seems to do better, I think
df <-as.data.frame(assay(rlog_dds)[,12:13])

ggplot(df, aes(x= RNAseq_OLD_rep12.bam, y=RNAseq_YOUNG_rep01.bam)) +geom_hex(bins=100, colour="orange", fill="black")+coord_fixed()+theme_classic()

#PCA
pca_plot<-plotPCA(rlog_dds)
ggsave("pca_RNAseq.pdf", pca_plot)

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
row.names(dds)[plot_load$var.names] 

```
