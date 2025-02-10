### Create Supplementary Table

### how does the data differ using the separate vs combined thresholds?

snps3<-read.csv("temp/genome/WG_CMHtest_results.csv")
indels3<-read.csv("temp/genome/indels_CMHtest_results.csv")
ann_snps<-read.table("data/genome/annotated_snps.txt", header=TRUE)
ann_indels<-read.table("data/genome/annotated_indels.txt", header= TRUE)

### REFORMAT annotations
## reformat ann to match working format
roman_to_chr <- c(
  "I" = "chr1", "II" = "chr2", "III" = "chr3", "IV" = "chr4", "V" = "chr5",
  "VI" = "chr6", "VII" = "chr7", "VIII" = "chr8", "IX" = "chr9", "X" = "chr10",
  "XI" = "chr11", "XII" = "chr12", "XIII" = "chr13", "XIV" = "chr14", "XV" = "chr15",
  "XVI" = "chr16", "Mito"= "chrmito"
)
# Replace Roman numerals in the CHROM column with chromosome notation
ann_snps <- ann_snps %>%
  mutate(CHROM = roman_to_chr[CHROM])

ann_indels <- ann_indels %>%
  mutate(CHROM = roman_to_chr[CHROM])

### Find the Combined threshold - use bonferroni correction for all data a=0.05
threshold_b_combined=0.05/(nrow(indels3)+nrow(snps3)) # this is a combined significance threshold
thresh_b_combined_log=-log10(threshold_b_combined)

### Applying combined threshold
combinedthresh_snps <- snps3[snps3$logp > thresh_b_combined_log, ] ## all snps above BF correction for combined SNPs/indels
combinedthresh_indels <- indels3[indels3$logp > thresh_b_combined_log, ] ## all indels above BF correction for combined SNPs/indels


SNP_NUC<-subset(combinedthresh_snps, CHROM!="chrmito") # nuclear SNPs
nrow(SNP_NUC) #624 nuclear SNPs
indel_nuc<-subset(combinedthresh_indels, CHROM!="chrmito")
nrow(indel_nuc) #104 nuclear indels

nrow(combinedthresh_snps) #669 SNPs
nrow(combinedthresh_indels) #132 indels

combinedthresh_snps$Type <- "SNP"
combinedthresh_indels$Type <- "Indel"

combined_sigs<-rbind(combinedthresh_snps, combinedthresh_indels)
nrow(combined_sigs) #801 combined SNPs/indels #includes mito

View(combined_sigs)

# Save the file to the new directory
write.csv(combined_sigs, file = "temp/genome/sigs_SNPs&indels_padj<0.05.csv", row.names = FALSE)

nrow(combined_sigs) ## 801 #sig list p<0.05 #no annotations #but need to merge annotation files with indels vs snps separately bc in different files
ann_snps<-merge(ann_snps, combinedthresh_snps, by=c("CHROM", "POS", "REF", "ALT"))
ann_indels<-merge(ann_indels, combinedthresh_indels, by=c("CHROM", "POS", "REF", "ALT"))

combined_ann_sigs<-rbind(ann_snps, ann_indels)
nrow(combined_ann_sigs) # missing 2 indels

missing_rows <- anti_join(combined_sigs, combined_ann_sigs, by = c("CHROM", "POS"))
View(missing_rows) # both in chr 3 #ARS 319 #Highly-active subtelomeric autonomously replicating sequence, initiates replication in ~90% of cell cycles, LOGP>17 ## located near RDS1 gene implicated in SNP data

View(combined_ann_sigs)
combined_ann_sigs2<-combined_ann_sigs[, -c(5:8, 13:61)]
View(combined_ann_sigs2)
ann_sigs_all<-combined_ann_sigs2[,c(1:8, 56:63)]

# Save the file to the new directory
write.csv(ann_sigs_all, file = "temp_tables/supp_table_complete_sig_list_p<0.05.csv", row.names = FALSE)

nrow(ann_sigs_all) #799 total sigs # missing 2 indels from chr3



