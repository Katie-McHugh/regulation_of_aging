###############################################################################
### Genic vs non-genic regions
###############################################################################

#------------------------------------------------------------------------------
# First run 00_common_gene_lists.R to generate datasets needed to load in

#------------------------------------------------------------------------------
### load in datasets

sigs_05_sep<-read.csv("temp/comparisons/sig_ALLannotations.csv")
sigs_05_f<-read.csv( "temp/comparisons/sig_FIRSTannotation.csv")

#------------------------------------------------------------------------------
#### What are the possible annotation categories? 

ann_cat<-unique(sigs_05_sep$Annotation)
print(ann_cat) #16 total categories...
## Genic categories: synonymous_variant, missense_variant, stop_gained, 
## start_lost, frameshift_variant, "frameshift_variant&start_lost", 
##"conservative_inframe_insertion", "disruptive_inframe_deletion", 
## "disruptive_inframe_insertion", "conservative_inframe_deletion" 
## "frameshift_variant&stop_lost&splice_region_variant"

## Non-genic: "intron_variant", "intergenic_region", "upstream_gene_variant",
## "downstream_gene_variant", "non_coding_transcript_exon_variant" 

#------------------------------------------------------------------------------
## remove any non-genic variants
cat_remove<- c("intron_variant", "intergenic_region", "upstream_gene_variant", 
               "downstream_gene_variant", "non_coding_transcript_exon_variant")
genic_sigs_sep<- sigs_05_sep[!sigs_05_sep$Annotation %in% cat_remove,]
nrow(genic_sigs_sep) #473 genic annotations

## find all intergenic variants
intergenic_sigs_f<- sigs_05_f[sigs_05_f$Annotation %in% cat_remove,]
nrow(intergenic_sigs_f) #336 intergenic annotations (df shows first annotation only)


intergenic_sigs_sep<- sigs_05_sep[sigs_05_sep$Annotation %in% cat_remove,]
View(intergenic_sigs_sep)
#View(intergenic_sigs_sep) to see all annotations in "Annotation" col

#------------------------------------------------------------------------------
### what are the intergenic variants? Various tables to summarize data

intergenic_annotations <- intergenic_sigs_f %>%
  group_by(Annotation) %>%
  summarize(Annotation_Count = n()) %>%
  ungroup() %>%
  mutate(Total_Annotations = sum(Annotation_Count)) %>%
  mutate(Percentage = (Annotation_Count / Total_Annotations) * 100) %>%
  ungroup()
View(intergenic_annotations)

intergenic_annotations_type <- intergenic_sigs_sep %>%
  group_by(Transcript_BioType) %>%
  summarize(Type_Count = n()) %>%
  ungroup() %>%
  mutate(Total_Annotations = sum(Type_Count)) %>%
  mutate(Percentage = (Type_Count / Total_Annotations) * 100) %>%
  ungroup()
View(intergenic_annotations_type)

transcript_biotype<-unique(intergenic_sigs_sep$Transcript_BioType)
print(transcript_biotype)
feature_type<-unique(intergenic_sigs_sep$Feature_Type)
print(feature_type)

# what about in the first annotation? 
genic_sigs_f<- sigs_05_f[!sigs_05_f$Annotation %in% cat_remove,]
nrow(genic_sigs_f) #463 genic annotations # very similar

nrow(intergenic_sigs_f) # PLUS 2 in ARS that weren't annotated

#463+336+2 # genic+ intergenic + 2 unannotated ARS

# most genic annotations are the first annotation, unsurprisingly
View(genic_sigs_f)
nrow(genic_sigs_sep)

#which ones differ? 
genic_diffs <- anti_join(genic_sigs_sep, genic_sigs_f, 
                         by = c("CHROM", "POS", "Annotation", "Gene_Name"))

nrow(genic_diffs) # this is the right number of rows
View(genic_diffs) ## the first annotation is still genic...there are just 
### multiple genic annotations (see note below)

genic_list<-unique(genic_sigs_f$Gene_Name)
length(genic_list)
print(genic_list)

## write list of genic SNPs for GO analysis
writeLines(genic_list, "temp/comparisons/genic_GO_list.txt")

#################### NOTES ON DIFFERENCES IN GENE LISTS ################

# Genic variants are typically the first annotation
## the ten "missing genes" are cases when there are TWO genic annotations for a 
## single gene

### Most likely, the first annotation is the more useful one, so keep that list

#------------------------------------------------------------------------------

write.csv(genic_sigs_f, file="temp/comparisons/genic_sigs.csv")
## more complete dataframe of significant genic SNPs

#------------------------------------------------------------------------------

## Create a combined GO-list for genic SNPs and RNA-seq genes
rna<-read.table("results/transcriptome/DGEgene_list.txt")
dna<-as.vector(genic_list)
rna<-as.vector(rna)
combined<-c(dna, rna)
combined<-as.character(combined)
writeLines(combined, "temp/comparisons/genic_combined_GO_list.txt")

#------------------------------------------------------------------------------

