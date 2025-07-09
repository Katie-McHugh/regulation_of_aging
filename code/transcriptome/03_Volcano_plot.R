#############################################################################
### Volcano Plot
#############################################################################

## Load Libraries
#------------------------------------------------------------------------------
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)
library(ggplot2)

#------------------------------------------------------------------------------

## Load in Data
#------------------------------------------------------------------------------
### Load normalized counts for visualization
res_adj<-read.csv("temp/transcriptome/rnaseq_results_batch_adjusted.csv", header=TRUE)
selected_genes_adj<-read.csv("temp/transcriptome/RNA_genes_p<0.1.csv", row.names= "X")
gene_key<- read.table("temp/transcriptome/key_geneIDtoName.txt")
View(res_adj)
## Simple Volcano Plot
#------------------------------------------------------------------------------
### Load list of genes for heatmap

res_adj2<-merge(res_adj, gene_key, by.x="X", by.y="Gene_ID")

selected_genes<-c("AHP1", "CCW12", "DIA4", "DNM1", "EKI1", "MRP1", "WSC4", 
                  "RFA3", "FIT2", "SPG5", "TOS2", "YDR543C", "YJL218W", "YNR066C")

labs <- ifelse(res_adj2$Gene_Name %in% selected_genes, res_adj2$Gene_Name, NA)


# ALL DEGs

custom_colors <- ifelse(res_adj2$log2FoldChange < 0, "steelblue", "red")
names(custom_colors) <- rownames(res_adj2)


p<-EnhancedVolcano(res_adj2,
                lab = labs,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'All DEGs',
                FCcutoff = FALSE,
                #col = c("green", "steelblue", "orange", "red"), #control point colors
                pCutoff = 0.1,
                colCustom = custom_colors,
                colAlpha = 0.8,
                cutoffLineCol = "grey45",
                pointSize = 3.0,
                labSize = 2.0, 
                selectLab = NA,
                xlim = c(-2, 2), 
                ylim=c(0, 3)
                )



plot(p)

selected_data <- subset(res_adj2, Gene_Name %in% selected_genes)

#------------------------------------------------------------------------------
### split data for easier labelling

# Split the selected data
neg_LFC <- selected_data[selected_data$log2FoldChange < 0, ]
pos_LFC <- selected_data[selected_data$log2FoldChange >= 0, ]
#------------------------------------------------------------------------------

# Add bigger hollow circles around those points on top of existing plot
p2<-p + 
  geom_point(data = selected_data, 
             aes(x = log2FoldChange, y = -log10(padj)),
             shape = 21, size = 3.5, color = "black", fill = NA, stroke = 0.5)+ 
  geom_text_repel(data = neg_LFC,
                  aes(x = log2FoldChange, y = -log10(padj), label = Gene_Name),
                  size = 4.1, 
                  fontface = "bold",
                  box.padding = 1, 
                  point.padding = 1,
                  segment.color = "black", 
                  segment.size = 0.2,
                  min.segment.length = 0.1, 
                  max.overlaps = Inf,
                  force = 2,
                  nudge_x = -0.1, 
                  nudge_y = 0.35) +
  geom_text_repel(data = pos_LFC,
                  aes(x = log2FoldChange, y = -log10(padj), label = Gene_Name),
                  size = 4.1, 
                  fontface = "bold",
                  box.padding = 1, 
                  point.padding = 1,
                  segment.color = "black", 
                  segment.size = 0.2,
                  min.segment.length = 0.1, 
                  max.overlaps = Inf,
                  force = 4,
                  nudge_x = 0.15, 
                  nudge_y = 0.35) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank())+
  labs(title = NULL, subtitle = NULL, caption = NULL) +  # Remove title, subtitle, caption
  theme(legend.position = "none")         


plot(p2)

#pdf("temp_figs/Volcano_DESEQadj_p<0.1.pdf", width = 8, height = 6)
pdf("figures/Figure3_Volcano.pdf", width = 8, height = 6)
p2
dev.off()

pdf("temp_figs/Volcano_DESEQp<0.1_v4.pdf", width = 8, height = 6)
p2
dev.off()


