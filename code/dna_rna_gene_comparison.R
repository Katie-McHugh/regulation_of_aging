### Analysis of implicated genomic annotations

## load in datasets

sigs_05_sep<-read.csv("temp/sig_ALLannotations.csv")
sigs_05_f<- read.csv("temp/sig_FIRSTannotation.csv")

# load in gene lists (genic vs non-genic script)
gene_list_all_ann<- read.table("temp/ALLannotations_list.txt", header=FALSE)
gene_list_first_ann<- read.table("temp/FIRSTannotation_list.txt", header=FALSE)
rna_list<-read.table("temp/DGEgene_list.txt", header= FALSE)
colnames(rna_list)[1]<-"Gene_Name"
View(gene_list_first_ann)

### Look at overlap between DNA and RNA lists

common_items <-intersect(rna_list$Gene_Name, gene_list_first_ann$V1)
print(common_items) ### WSC4, RFA3

#### Restrictive list for RNA, ALL annotations for DNA
common_items2 <-intersect(rna_list$Gene_Name, gene_list_all_ann$V1)
print(common_items2) #"YDR543C" "EKI1"    "YNR066C" "WSC4"    "FIT2"    "RFA3"    "DNM1"    "YJL218W"
### EKI1, DNM1 are downstream gene variants

### View common items

### This shows just the first annotation (only RFA3 and wsc4)
sigs_05_both<-sigs_05_f[sigs_05_f$Gene_Name %in% common_items,]
View(sigs_05_both)

### This shows all annotations
sigs_05_common<-sigs_05_sep[sigs_05_sep$Gene_Name %in% common_items2,]
View(sigs_05_common)

write.csv(sigs_05_both, file="temp/shared_genes_FIRSTann.csv")
write.csv(sigs_05_common, file="temp/shared_genes_ALLann.csv")

## Now look at variant identities
gene_annotation_table <- sigs_05_common %>%
  group_by(Gene_Name) %>%
  summarize(Annotations = paste(sort(unique(Annotation)), collapse = ", ")) %>%
  ungroup()
View(gene_annotation_table)

# Create a new table that counts how many genes have each combination of annotations
gene_annotation_count <- sigs_05_common %>%
  group_by(Gene_Name, Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup()
View(gene_annotation_count)

gene_annotation_count_FIRST <- sigs_05_both %>%
  group_by(Gene_Name, Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup()
View(gene_annotation_count_FIRST)

