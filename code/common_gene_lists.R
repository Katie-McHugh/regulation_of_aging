### Comparing genomic and transcriptomic data to find shared genes

# load in significant genome and transcriptome data
sigs_05<-read.csv("temp_tables/supp_table_complete_sig_list_p<0.05.csv")
rna_list<-read.csv("temp/DGEgene_list.csv")

### Look at the annotations for our gene list

#### Start by just looking at the first (and presumably most relevant) annotation for each SNP/indel

##### Step 1: Extract the first annotation by splitting at the first comma
sigs_05_f <- sigs_05 %>%
  mutate(First_ANN = sapply(strsplit(ANN, ","), `[`, 1))

##### Step 2: Separate the First_ANN column by pipe character into multiple columns
sigs_05_f <- sigs_05_f %>%
  separate(First_ANN, into = c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "INFO"), sep = "\\|", fill = "right")

nrow(sigs_05_f) #799 SNPs/indels total
View(sigs_05_f)

### If we wanted to look at ALL annotations (not just the first one in SNPeff)
#### Step 1: Split the ANN column into individual annotations ### this creates a new row for each annotation
sigs_05_all <- sigs_05 %>%
  mutate(Annotations_List = strsplit(as.character(ANN), ",")) %>%  # Split the annotations by comma
  unnest(Annotations_List)  # Expand the list into multiple rows

#### Step 2: Separate the components of each annotation
sigs_05_sep <- sigs_05_all %>%
  separate(Annotations_List, 
           into = c("Indel", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", 
                    "Feature_Type", "Feature_ID", "Transcript_BioType", "Coding_Change", 
                    "Amino_Acid_Change"), 
           sep = "\\|", fill = "right", extra = "merge")

nrow(sigs_05_sep) #5510 total annotations for the 799 SNP/indel sites

gene_list_all_ann<-unique(sigs_05_sep$Gene_Name) #1180 genes within ALL annotations
length(gene_list_all_ann)

gene_list_first_ann<-unique(sigs_05_f$Gene_Name) #genes within the FIRST annotation
length(gene_list_first_ann) #233 features ### this includes non-protein coding regions

### Create text file with lists for GO-analysis and comparisons
writeLines(gene_list_all_ann, "temp/ALLannotations_list.txt")
writeLines(gene_list_first_ann, "temp/FIRSTannotation_list.txt")

### compare RNA and DNA lists
gene_list_all_ann<-as.data.frame(gene_list_all_ann)
gene_list_first_ann<-as.data.frame(gene_list_first_ann)

common_items <-intersect(rna_list$Gene_Name, gene_list_first_ann$gene_list_first_ann)
print(common_items) ### WSC4, RFA3

#### Restrictive list for RNA, ALL annotations for DNA
common_items2 <-intersect(rna_list$Gene_Name, gene_list_all_ann$gene_list_all_ann)
print(common_items2) #"YDR543C" "EKI1"    "YNR066C" "WSC4"    "FIT2"    "RFA3"    "DNM1"    "YJL218W"
### EKI1, DNM1 are downstream gene variants

### View comon items

### This shows just the first annotation (only RFA3 and wsc4)
sigs_05_both<-sigs_05_f[sigs_05_f$Gene_Name %in% common_items,]
View(sigs_05_both)

### This shows all annotations
sigs_05_common<-sigs_05_sep[sigs_05_sep$Gene_Name %in% common_items2,]
View(sigs_05_common)

write.csv(sigs_05_both, file="temp/shared_genes_FIRSTann.csv")
write.csv(sigs_05_common, file="temp/shared_genes_ALLann.csv")
