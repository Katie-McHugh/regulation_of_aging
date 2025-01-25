### Genic vs non-genic regions

### load in data
sigs_05_sep<-read.csv("temp/sig_ALLannotations.csv")
sigs_05_f<-read.csv( "temp/sig_FIRSTannotation.csv")

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

## remove any non-genic variants
cat_remove<- c("intron_variant", "intergenic_region", "upstream_gene_variant", "downstream_gene_variant", "non_coding_transcript_exon_variant")
genic_sigs_sep<- sigs_05_sep[!sigs_05_sep$Annotation %in% cat_remove,]
nrow(genic_sigs_sep) #473 genic annotations
View(genic_sigs_sep)
# what about in the first annotation? 
genic_sigs_f<- sigs_05_f[!sigs_05_f$Annotation %in% cat_remove,]
nrow(genic_sigs_f) #463 genic annotations # very similar
# most genic annotations are the first annotation, unsurprisingly

#which ones differ? 
genic_diffs <- anti_join(genic_sigs_sep, genic_sigs_f, by = c("CHROM", "POS", "Annotation", "Gene_Name"))
nrow(genic_diffs) # this is the right number of rows
View(genic_diffs) ## these all seem legit, so keep the longer list

write.csv(genic_sigs_sep, file="temp/genic_sigs.csv")
