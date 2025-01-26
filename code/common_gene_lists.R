### Creating simplified gene lists for next steps 
 ##### next steps = dna_rna_gene_comparison.R and genic_vs_nongenic.R scripts

# load in significant genome and transcriptome data
sigs_05<-read.csv("temp_tables/supp_table_complete_sig_list_p<0.05.csv")

### Look at the annotations for our DNA gene lists

#### Start by just looking at the first (and presumably most relevant) annotation for each SNP/indel

##### Step 1: Extract the first annotation by splitting at the first comma
sigs_05_f <- sigs_05 %>%
  mutate(First_ANN = sapply(strsplit(ANN, ","), `[`, 1))

##### Step 2: Separate the First_ANN column by pipe character into multiple columns
sigs_05_f <- sigs_05_f %>%
  separate(First_ANN, into = c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "INFO"), sep = "\\|", fill = "right")

nrow(sigs_05_f) #799 SNPs/indels total #matches length of original list

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
gene_list_first_ann<-unique(sigs_05_f$Gene_Name) #features within the FIRST annotation
length(gene_list_first_ann) #233 features ### this includes non-protein coding regions

### Create text file with lists for GO-analysis and comparisons #this is still only genomic data
write.csv(sigs_05_sep, "temp/sig_ALLannotations.csv")
write.csv(sigs_05_f, "temp/sig_FIRSTannotation.csv")

writeLines(gene_list_all_ann, "temp/ALLannotations_list.txt")
writeLines(gene_list_first_ann, "temp/FIRSTannotation_list.txt")

######

## gene list without mitochondria

dna<-subset(sigs_05_f, CHROM!="chrmito")
dna_list<-unique(dna$Gene_Name)
View(as.data.frame(dna_list))
rna<-unique(rna_list$Gene_Name)
combined<-unique(c(dna_list, rna))
View(combined)

## go doesn't like these gene names, wants the SGID instead
lookup_table <- data.frame(
  Gene_Name = c("KIN3", "KIN4", "RDS1", "BUL1", "NSP1", "TFS1", "TIP1", "TIP20", "HAP2", "IKI1"),
  SGDID = c("S000000071", "S000005759", "S000000703", "S000004888", 
            "S000003577", "S000004168", "S000000271", "S000003113", 
            "S000003206", "S000001230")
)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# Convert combined list to a data frame
combined_df <- data.frame(Gene_Name = combined)

# Perform a left join with the lookup table
result <- combined_df %>%
  left_join(lookup_table, by = "Gene_Name") %>%
  mutate(Final_Name = ifelse(is.na(SGDID), Gene_Name, SGDID)) %>%
  select(Final_Name) %>%
  pull()

writeLines(result, "temp/GO_combined_no_mt.txt")



