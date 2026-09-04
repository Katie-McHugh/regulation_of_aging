### Comparing Global Expression of DEGs vs non-DEGs
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

### Load results object
res_adj<-read.csv("data/clean/rnaseq_results_batch_adjusted_young_ref.csv", header=TRUE)

## just read in table instead of generating it

gene_key<- read.table("temp/transcriptome/key_geneIDtoName.txt")

head(res_adj)
sig_genes_adj_1<-merge(res_adj, gene_key, by.x="X", by.y="Gene_ID", all.x=TRUE)
sig_genes_adj_1$Gene_Name[res_adj$X == "YLR050C"] <- "EMA19" ## Name recognized in SGD

res_adj<-merge(res_adj, gene_key, by.x="X", by.y="Gene_ID", all.x=TRUE)

res_adj$Gene_Name[res_adj$X == "YLR050C"] <- "EMA19" ## Name recognized in SGD
res_adj$Gene_Name[res_adj$X == "YCR015C"] <- "CTO1" ## Name recognized in SGD

## rename columns
colnames(res_adj)[colnames(res_adj) == "X"] <- "Gene_ID"
#colnames(sig_genes_adj_1)[colnames(res_1) == "X"] <- "Gene_ID"

#------------------------------------------------------------------------------
## write to files

write.csv(res_adj, file="results/transcriptome/RNA_genes_all.csv") 

DEGs<- res_adj[which(res_adj$padj <= 0.1), ] 
non_DEGs<-res_adj[which(res_adj$padj >= 0.1), ] 

DEGs$LBM<-log10(DEGs$baseMean)
non_DEGs$LBM<-log10(non_DEGs$baseMean)

library(ggplot2)
library(patchwork)

p1 <- ggplot(DEGs, aes(x = LBM)) +
  geom_histogram(bins = 20, fill = "darkslategray4", alpha =0.7) +
  labs(x = "log10(baseMean)", y = "DEGs")

p2 <- ggplot(non_DEGs, aes(x = LBM)) +
  geom_histogram(bins = 20, fill = "chocolate3", alpha = 0.7) +
  labs(x = "log10(baseMean)", y = "non-DEGs")

hist_patch <- p1 / p2 +
  plot_layout(axis_titles = "collect") &
  scale_x_continuous(limits = range(c(DEGs$LBM, non_DEGs$LBM))) &
                       theme(
                         axis.title = element_text(size = 20),
                         axis.text = element_text(size = 16))
hist_patch

ggsave(filename = "figures/global_expression_comparison_2.png", plot = hist_patch, width = 7, height = 8, dpi = 300)

## alternate

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


mean(non_DEGs$LBM)
mean(DEGs$LBM)

## Wilcoxon test? 

?wilcox.test
head(DEGs)
head(non_DEGs)

res<-wilcox.test(DEGs$baseMean, non_DEGs$baseMean,
            alternative = "two.sided")
res
z_value <- qnorm(res$p.value / 2)
print(abs(z_value))
z_effect<- z_value/sqrt(nrow(DEGs)+nrow(non_DEGs))
z_effect

install.packages("rcompanion")
library(rcompanion)

# Works directly on vector variables
wilcoxonZ(DEGs$baseMean, non_DEGs$baseMean, paired = FALSE)

## technically there is a significant difference (p=0.004) 
## but the effect size is negligible (0.03)