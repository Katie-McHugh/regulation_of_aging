### Use DESeq object to find DEGs

### Load results object
res_adj<-read.csv("temp/rnaseq_results_batch_adjusted.csv", header=TRUE)

# find genes below p<0.05
sig_genes_adj <- res_adj[which(res_adj$padj <= 0.05), ] #27 genes 

#p<0.1
sig_genes_permissive_adj <- res_adj[which(res_adj$padj <= 0.1), ] #60 genes #includes RFA3 and FIT2...
nrow(sig_genes_permissive_adj)

ann<-read.table("data/annotated_snps.txt", header = TRUE)
roman_to_chr <- c(
  "I" = "chr1", "II" = "chr2", "III" = "chr3", "IV" = "chr4", "V" = "chr5",
  "VI" = "chr6", "VII" = "chr7", "VIII" = "chr8", "IX" = "chr9", "X" = "chr10",
  "XI" = "chr11", "XII" = "chr12", "XIII" = "chr13", "XIV" = "chr14", "XV" = "chr15",
  "XVI" = "chr16", "Mito"= "chrmito"
)
# Replace Roman numerals in the CHROM column with chromosome notation
ann <- ann %>%
  mutate(CHROM = roman_to_chr[CHROM])

ann<-as.data.frame(ann)
View(an)
ann_all <- ann %>%
  mutate(Annotations_List = strsplit(as.character(ann), ",")) %>%  # Split the annotations by comma
  unnest(Annotations_List)  # Expand the list into multiple rows

#### Step 2: Separate the components of each annotation
ann_sep <- ann_all %>%
  separate(Annotations_List, 
           into = c("Indel", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", 
                    "Feature_Type", "Feature_ID", "Transcript_BioType", "Coding_Change", 
                    "Amino_Acid_Change"), 
           sep = "\\|", fill = "right", extra = "merge")

## save separated out annotations file for later
write.table(ann_sep, file="temp/annotations_snps_sep.txt")

### create key for gene ID to gene name conversion
gene_key<-ann_sep[,5:6]
gene_key<-unique(gene_key)

## save key to temp folder
write.table(gene_key, file="temp/key_geneIDtoName.txt")

### Create list of significant genes

sig_genes_adj_1<-merge(sig_genes_adj, gene_key, by.x="X", by.y="Gene_ID", all.x=TRUE)
View(sig_genes_adj_1)
sig_genes_adj_1$Gene_Name[sig_genes_adj_1$X == "YLR050C"] <- "EMA19" ## Name recognized in SGD

sig_genes_permissive_adj_1<-merge(sig_genes_permissive_adj, gene_key, by.x="X", by.y="Gene_ID", all.x=TRUE)
View(sig_genes_permissive_adj_1)
sig_genes_permissive_adj_1$Gene_Name[sig_genes_permissive_adj_1$X == "YLR050C"] <- "EMA19" ## Name recognized in SGD
sig_genes_permissive_adj_1$Gene_Name[sig_genes_permissive_adj_1$X == "YCR015C"] <- "CTO1" ## Name recognized in SGD

## rename columns
colnames(sig_genes_permissive_adj_1)[colnames(sig_genes_permissive_adj_1) == "X"] <- "Gene_ID"
colnames(sig_genes_adj_1)[colnames(sig_genes_adj_1) == "X"] <- "Gene_ID"

write.csv(sig_genes_adj_1, file="temp/RNA_genes_p<0.05.csv") 
write.csv(sig_genes_permissive_adj_1, file="temp/RNA_genes_p<0.1.csv") 

gene_names_p05<-unique(sig_genes_adj_1$Gene_Name)
gene_names_p10<-unique(sig_genes_permissive_adj_1$Gene_Name)

writeLines(gene_names_p10, "temp/DGEgene_list.txt")


