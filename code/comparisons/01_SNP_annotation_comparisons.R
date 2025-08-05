### Comparing annotations for significant genes vs the whole genome
### SNPs only

#-------------------------------------------------------------------------------
#### Script reorganized annotated SNP data to look for evidence of overrepresented 
## categories in significant SNPs - 
##output us used for Fishers Test in 02_Fisherstest_annotations.R
#-------------------------------------------------------------------------------

## Load in data

### significant data
sigs_05_f<-read.csv("temp/comparisons/sig_FIRSTannotation.csv") 
sigs_05_f<-subset(sigs_05_f, TYPE=="SNP")
sigs_05_sep<-read.csv("temp/comparisons/sig_ALLannotations.csv")
sigs_05_sep<-subset(sigs_05_sep, TYPE=="SNP")

### whole genome annotations
ann_s<-read.table("data/raw/annotated_snps.txt", header=TRUE)  #snps only
#-------------------------------------------------------------------------------

###############################################################################
### Now look at SNPeffs whole-genome annotations
##############################################################################

################ edit this to look at SNPs and indels ########################

### Just look at SNPs - can change to be both
##ann<-rbind(ann_s, ann_i)
ann<-ann_s

#### Isolate first annotation
ann_f <- ann %>%
  mutate(First_ANN = sapply(strsplit(ANN, ","), `[`, 1))

#### Step 2: Separate the First_ANN column by pipe character into multiple columns
ann_f <- ann_f %>%
  separate(First_ANN, into = c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "INFO"), sep = "\\|", fill = "right")

### subset relevant information
ann_f2<-ann_f[,-c(13:60)]
ann_f3<-ann_f[c("CHROM", "POS", "TYPE", "Gene_ID", "Gene_Name", "Annotation","Annotation_Impact","Feature_Type", "Transcript_BioType", "EVENTLENGTH", "ANN")]


### first annotation for every annotated variant in the genome
gene_annotation_count_first_genome <- ann_f3 %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup()

### same but convert to percentages instead of counts
ref_first_ann <- ann_f3 %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup() %>%
  mutate(Total_Annotations = sum(Annotation_Count)) %>%
  mutate(Percentage = (Annotation_Count / Total_Annotations) * 100) %>%
  ungroup()

write.csv(ref_first_ann, "temp/comparisons/ref_annotations_count.csv")

#-------------------------------------------------------------------------------

### Looking at first annotation for all sigs...what is the breakdown?
gene_annotation_count <- sigs_05_f %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup()

write.csv(gene_annotation_count, file = "temp/comparisons/significant_annotations_count.csv", row.names = FALSE, quote=FALSE) 



