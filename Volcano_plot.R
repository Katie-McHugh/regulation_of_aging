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
View(norm_dds)
norm_dds<-read.csv("temp/transcriptome/normalized_counts_deseq.csv")
res_adj<-read.csv("temp/transcriptome/rnaseq_results_batch_adjusted.csv", header=TRUE)
selected_genes_adj<-read.csv("temp/transcriptome/RNA_genes_p<0.1.csv", row.names= "X")
gene_key<- read.table("temp/transcriptome/key_geneIDtoName.txt")

## Simple Volcano Plot
#------------------------------------------------------------------------------
### Load list of genes for headmap

res_adj2<-merge(res_adj, gene_key, by.x="X", by.y="Gene_ID")

selected_genes<-c("AHP1", "CCW12", "DIA4", "DNM1", "EKI1", "MRP1", "WSC4", 
                  "RFA3", "FIT2", "SPG5", "TOS2", "YDR543C", "YJL218W", "YNR066C")

labs <- ifelse(res_adj2$Gene_Name %in% selected_genes, res_adj2$Gene_Name, NA)

# ALL DEGs

p<-EnhancedVolcano(res_adj2,
                lab = labs,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'All DEGs',
                FCcutoff = TRUE,
                #col = c("red", "grey70", "red", "red"), #control point colors
                pCutoff = 0.1,
                cutoffLineCol = "grey50",
                pointSize = 3.0,
                labSize = 2.0, 
                selectLab = NA,
                xlim = c(-2, 2), 
                ylim=c(0, 3)
                )

selected_data <- subset(res_adj2, Gene_Name %in% selected_genes)

# Add bigger hollow circles around those points on top of existing plot
p2<-p + 
  geom_point(data = selected_data, 
             aes(x = log2FoldChange, y = -log10(padj)),
             shape = 21, size = 3.5, color = "black", fill = NA, stroke = 0.5)+ 
  geom_text_repel(data = selected_data,
                  aes(x = log2FoldChange, y = -log10(padj), label = Gene_Name),
                  size = 3,
                  fontface = "bold",
                  box.padding = 0.8,          # More space around label box
                  point.padding = 1,        # More space around data point
                  segment.color = "black",    # Keep connectors
                  segment.size = 0.3,         # Thinner connector lines
                  min.segment.length = 0.2,   # Avoid tiny/tangled segments
                  max.overlaps = Inf,         # Avoid skipping labels due to overlap
                  force = 3,                # Stronger repulsion (try 1–5)
                  nudge_y = 0.1              # Optional vertical nudge 
  ) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank())+
  labs(title = NULL, subtitle = NULL, caption = NULL) +  # Remove title, subtitle, caption
  theme(legend.position = "none")         


plot(p2)

#pdf("temp_figs/Volcano_DESEQadj_p<0.1.pdf", width = 8, height = 6)
pdf("temp_figs/Volcano_DESEQp<0.1.pdf", width = 4, height = 3)
p2
dev.off()

pdf("temp_figs/Volcano_DESEQp<0.1_v2.pdf", width = 4, height = 3)
p2
dev.off()

# P<0.1
EnhancedVolcano(selected_genes_adj,
                lab = selected_genes_adj$Gene_Name,
                x = 'log2FoldChange',
                y = 'padj')


library(ggplot2)

# Base volcano points
ggplot(res_adj2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(color = "grey60", shape = 20) +  # all points
  
  # Add highlighted points with different shape and color
  geom_point(data = subset(res_adj2, Gene_Name %in% selected_genes),
             aes(x = log2FoldChange, y = -log10(padj)),
             color = "red", shape = 21, size = 3) +
  
  # Add labels for selected genes
  geom_text(data = subset(res_adj2, Gene_Name %in% selected_genes),
            aes(label = Gene_Name), size = 3, vjust = 1.5)

