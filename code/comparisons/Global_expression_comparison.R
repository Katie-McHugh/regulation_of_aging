### Comparing Global Expression of DEGs vs non-DEGs
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

### Load results object
res_adj<-read.csv("data/clean/rnaseq_results_batch_adjusted_young_ref.csv", header=TRUE)

## just read in table instead of generating it

gene_key<- read.table("temp/transcriptome/key_geneIDtoName.txt")


sig_genes_adj_1<-merge(sig_genes_adj, gene_key, by.x="X", by.y="Gene_ID", all.x=TRUE)
sig_genes_adj_1$Gene_Name[res_adj$X == "YLR050C"] <- "EMA19" ## Name recognized in SGD

res_adj<-merge(res_adj, gene_key, by.x="X", by.y="Gene_ID", all.x=TRUE)

res_adj$Gene_Name[res_adj$X == "YLR050C"] <- "EMA19" ## Name recognized in SGD
res_adj$Gene_Name[res_adj$X == "YCR015C"] <- "CTO1" ## Name recognized in SGD

## rename columns
colnames(res_adj)[colnames(res_adj) == "X"] <- "Gene_ID"
colnames(sig_genes_adj_1)[colnames(res_1) == "X"] <- "Gene_ID"

#------------------------------------------------------------------------------
## write to files

write.csv(res_adj, file="results/transcriptome/RNA_genes_all.csv") 

DEGs <- res_adj[which(res_adj$padj <= 0.1), ] 
non_DEGs<-res_adj[which(res_adj$padj >= 0.1), ] 

DEGs$LBM<-log10(DEGs$baseMean)
non_DEGs$LBM<-log10(non_DEGs$baseMean)

library(ggplot2)
library(patchwork)

p1 <- ggplot(DEGs, aes(x = LBM)) +
  geom_histogram(bins = 10, fill = "darkslategray4", alpha =0.7) +
  labs(x = "log10(baseMean)", y = "DEGs")

p2 <- ggplot(non_DEGs, aes(x = LBM)) +
  geom_histogram(bins = 10, fill = "chocolate3", alpha = 0.7) +
  labs(x = "log10(baseMean)", y = "non-DEGs")

hist_patch <- p1 / p2 +
  plot_layout(axis_titles = "collect") &
  scale_x_continuous(limits = range(c(DEGs$LBM, non_DEGs$LBM))) &
                       theme(
                         axis.title = element_text(size = 20),
                         axis.text = element_text(size = 16))
hist_patch

ggsave(filename = "figures/global_expression_comparison.png", plot = hist_patch, width = 7, height = 8, dpi = 300)


View(DEGs)

?hist

DEG_hist<-hist(DEGs$LBM, 
     breaks = 20, 
     xlab = "log10(baseMean)", 
     col = "darkslategray4",
     main = NA)

non_DEG_hist<-hist(non_DEGs$LBM, 
     breaks = 20, 
     xlab = "log10(baseMean)", 
     col = "chocolate3",
     main = NA)

library(patchwork)

hist_patch<- ( DEG_hist / non_DEG_hist)

mean(non_DEGs$baseMean)
mean(DEGs$baseMean, 
     )
