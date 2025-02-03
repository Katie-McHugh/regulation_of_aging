### Comparing annotations for significant genes vs the whole genome
### SNPs only

######################load in data############################################

### significant data
sigs_05_f<-read.csv("temp/sig_FIRSTannotation.csv")
sigs_05_f<-subset(sigs_05_f, TYPE=="SNP")
sigs_05_sep<-read.csv("temp/sig_ALLannotations.csv")
sigs_05_sep<-subset(sigs_05_sep, TYPE=="SNP")
### whole genome annotations
###ann_i<-read.table("data/annotated_indels.txt", header=TRUE)
ann_s<-read.table("data/annotated_snps.txt", header=TRUE)

#############################################################################

### Create tables to visualize the data

### first annotation in SNPeff for given variant, sorted by annotation type
gene_annotation_table_first<- sigs_05_f %>%
  group_by(Gene_Name) %>%
  summarize(Annotations = paste(sort(unique(Annotation)), collapse = ", ")) %>%
  ungroup()


### Create a new table that counts how many genes have each combination of annotations
gene_annotation_count_first <- sigs_05_f %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup()

###########################################
# use gene_annotation_count_first for stats
############################################


### Convert counts to percentages
gene_annotation_first_percent <- sigs_05_f %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup() %>%
  mutate(Total_Annotations = sum(Annotation_Count)) %>%
  mutate(Percentage = (Annotation_Count / Total_Annotations) * 100) %>%
  ungroup()


### Looking at all annotations from significant data 
### _sep = separated out each variants annotations into unique rows

gene_annotation_count_all <- sigs_05_sep %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup()

### Percentages
gene_annotation_all_percent <- sigs_05_sep %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup() %>%
  mutate(Total_Annotations = sum(Annotation_Count)) %>%
  mutate(Percentage = (Annotation_Count / Total_Annotations) * 100) %>%
  ungroup()


###############################################################################
### Now look at SNPeffs whole-genome annotations
##############################################################################

################ edit this to look at SNPs and indels ########################

### Just look at SNPs - can change to be both
##ann<-rbind(ann_s, ann_i)
ann<-ann_s

#############################################################################


#### Step 1: Extract the first annotation by splitting at the first comma
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

write.csv(ref_first_ann, "ref_annotations_count.csv")

View(ref_first_ann)

###########################################
# use ref_first_ann for stats
############################################


### this looks at ALL annotations, not just the first annotation

#### Step 1: Split the ANN column into individual annotations ### this creates a new row for each annotation
ann_all <- ann %>%
  mutate(Annotations_List = strsplit(as.character(ANN), ",")) %>%  # Split the annotations by comma
  unnest(Annotations_List)  # Expand the list into multiple rows

#### Step 2: Separate the components of each annotation
ann_sep <- ann_all %>%
  separate(Annotations_List, 
           into = c("Indel", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", 
                    "Feature_Type", "Feature_ID", "Transcript_BioType", "Coding_Change", 
                    "Amino_Acid_Change"), 
           sep = "\\|", fill = "right", extra = "merge")

ann_sep<-ann_sep[c("CHROM", "POS", "TYPE", "Gene_ID", "Gene_Name", "Annotation","Annotation_Impact","Feature_Type", "Transcript_BioType", "EVENTLENGTH", "ANN")]


### all annotations for all annotated variants
gene_annotation_all_genome <- ann_sep %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup()

##Look at percentages
ref_percent_all <- ann_sep %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup() %>%
  mutate(Total_Annotations = sum(Annotation_Count)) %>%
  mutate(Percentage = (Annotation_Count / Total_Annotations) * 100) %>%
  ungroup()

# Select only Annotation and Percentage columns, rename the Percentage column 
#to clarify which is which
gene_annotation_first_percent<-as.data.frame(gene_annotation_first_percent)
colnames(gene_annotation_first_percent)
class(ref_first_ann)

table_1 <- gene_annotation_first_percent[c("Annotation", "Percentage")]
colnames(table_1)[colnames(table_1) == "Percentage"] <- "Data_first"

table_2 <- ref_first_ann[c("Annotation", "Percentage")]
colnames(table_2)[colnames(table_2) == "Percentage"] <- "Ref_first"

table_3 <- gene_annotation_all_percent[c("Annotation", "Percentage")]
colnames(table_3)[colnames(table_3) == "Percentage"] <- "Data_all"

table_4 <- ref_percent_all[c("Annotation", "Percentage")]
colnames(table_4)[colnames(table_4) == "Percentage"] <- "Ref_all"


# Perform full join to combine all tables by the Annotation column
combined_table <- table_1 %>%
  full_join(table_2, by = "Annotation") %>%
  full_join(table_3, by = "Annotation") %>%
  full_join(table_4, by = "Annotation")

# Replace NA with "0%" in all scenario columns
combined_table <- combined_table %>%
  mutate(across(starts_with("Scenario"), ~ ifelse(is.na(.), 0, .)))

combined_table <- combined_table %>%
  mutate(across(where(is.numeric), ~ formatC(., format = "f", digits = 5)))

View(combined_table)

write.csv(combined_table, file = "temp/annotation_comparisons.csv", row.names = FALSE, quote=FALSE) 


### table with just SNPs? 
## combine annotations to just compare genic vs non-genic? or do that as well> 

###################################################################################

### Looking at first annotation for all sigs...what is the breakdown?
gene_annotation_count <- sigs_05_f %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup()

write.csv(gene_annotation_count, file = "temp/significant_annotations_count.csv", row.names = FALSE, quote=FALSE) 

### also make one just with SNP data
View(sigs_05_f)
sigs_05_f_snps<-subset(sigs_05_f, TYPE=="SNP")

gene_annotation_count_snps <- sigs_05_f_snps %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup()
View(gene_annotation_count_snps)

write.csv(gene_annotation_count, file = "temp/significant_annotations_count.csv", row.names = FALSE, quote=FALSE) 
write.csv(gene_annotation_count_snps, file = "temp/significant_annotations_count_snps_only.csv", row.names = FALSE, quote=FALSE) 

# Convert Data_first to numeric and replace "NA" strings with actual NA
combined_table <- combined_table %>%
  mutate(Data_first = na_if(Data_first, "NA"),  # Replace "NA" string with actual NA
         Data_first = as.numeric(Data_first))  # Convert to numeric

# Now check for NA values and extract annotations
na_categories <- combined_table %>% 
  filter(is.na(Data_first)) %>% 
  pull(Annotation) %>% 
  unique()  # Get unique annotations with NA

# View the result
print(na_categories)

# Combine specified categories into "Other"
combined_table2 <- combined_table %>%
  mutate(Annotation = if_else(Annotation %in% na_categories, "Other", Annotation))

# Group by Annotation and summarize numeric columns
combined_table2 <- combined_table2 %>%
  group_by(Annotation) %>%
  summarize(
    Data_first = sum(as.numeric(Data_first), na.rm = TRUE),
    Ref_first = sum(as.numeric(Ref_first), na.rm = TRUE),
    Data_all = sum(as.numeric(Data_all), na.rm = TRUE),
    Ref_all = sum(as.numeric(Ref_all), na.rm = TRUE),
    .groups = 'drop'  # Ungroup the result
  ) %>%
  arrange(factor(Annotation, levels = c(setdiff(unique(Annotation), "Other"), "Other")))

View(combined_table2)

write.csv(combined_table2, file = "temp/Percent_annotations_comparison.csv", row.names = FALSE, quote=FALSE) 
