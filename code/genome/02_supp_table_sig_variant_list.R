### Create Supplementary Table

library(dplyr)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
## read in data
snps3<-read.csv("results/genome/WG_CMHtest_results.csv")
indels3<-read.csv("results/genome/indels_CMHtest_results.csv")
ann_snps<-read.table("data/raw/annotated_snps.txt", header=TRUE)
ann_indels<-read.table("data/raw/annotated_indels.txt", header= TRUE)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

### Find the Combined threshold - use bonferroni correction for all data a=0.05

# this is a combined significance threshold
threshold_b_combined=0.05/(nrow(indels3)+nrow(snps3)) 
thresh_b_combined_log=-log10(threshold_b_combined)

### Applying combined threshold

## all snps above BF correction for combined SNPs/indels
combinedthresh_snps <- snps3[snps3$logp > thresh_b_combined_log, ] 

## all indels above BF correction for combined SNPs/indels
combinedthresh_indels <- indels3[indels3$logp > thresh_b_combined_log, ] 

#------------------------------------------------------------------------------

# nuclear SNPs
SNP_NUC<-subset(combinedthresh_snps, CHROM!="chrmito") 
nrow(SNP_NUC) #624 nuclear SNPs

## nuclear indels
indel_nuc<-subset(combinedthresh_indels, CHROM!="chrmito")
nrow(indel_nuc) #104 nuclear indels

nrow(combinedthresh_snps) #669 SNPs
nrow(combinedthresh_indels) #132 indels

combinedthresh_snps$Type <- "SNP"
combinedthresh_indels$Type <- "Indel"

combined_sigs<-rbind(combinedthresh_snps, combinedthresh_indels)
nrow(combined_sigs) #801 combined SNPs/indels #includes mito

#------------------------------------------------------------------------------

# Save the file to the new directory
write.csv(combined_sigs, file = "results/genome/sigs_SNPs&indels_padj<0.05.csv", row.names = FALSE)

#------------------------------------------------------------------------------

nrow(combined_sigs) ## 801 #sig list p<0.05 
#no annotations #but need to merge annotation files with indels vs snps separately bc in different files
ann_snps<-merge(ann_snps, combinedthresh_snps, by=c("CHROM", "POS", "REF", "ALT"))
ann_indels<-merge(ann_indels, combinedthresh_indels, by=c("CHROM", "POS", "REF", "ALT"))

combined_ann_sigs<-rbind(ann_snps, ann_indels)
nrow(combined_ann_sigs) # missing 2 indels

missing_rows <- anti_join(combined_sigs, combined_ann_sigs, by = c("CHROM", "POS"))
View(missing_rows) # both in chr 3 
#ARS 319 #Highly-active subtelomeric autonomously replicating sequence, 
#initiates replication in ~90% of cell cycles, LOGP>17 
## located near RDS1 gene implicated in SNP data

#------------------------------------------------------------------------------
# Write as file
write.csv(missing_rows, file="temp/genome/unannotated_indels.csv", row.names = FALSE)
#------------------------------------------------------------------------------

combined_ann_sigs2<-combined_ann_sigs[, -c(5:8, 13:61)]
ann_sigs_all<-combined_ann_sigs2[,c(1:8, 56:63)]
ann_sigs_all<-rbind(ann_sigs_all)

# Save the file to the new directory
write.csv(ann_sigs_all, file = "results/genome/supp_table_complete_sig_list_p<0.05.csv", row.names = FALSE)

nrow(ann_sigs_all) #799 total sigs # missing 2 indels from chr3

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
## Merge summary table and SNP table


# SNP table with allele counts
head(combined_sigs)

# reorder col names
new_order <- c("CHROM", "POS", "REF", "ALT", "TYPE", "ANN", "logp", "pval", "mincov", "maxcov", "sum_mac", "sum_cov", "maf", "EVENTLENGTH", "MULTI.ALLELIC", "Type")
ann_sigs_new <- ann_sigs_all[ , new_order ]
ann_sigs_new<-ann_sigs_new[, 1:13]


# Summary table with annotations
head(ann_sigs_all)

## merge together
merged_sig_table <-merge(ann_sigs_new, combined_sigs, by.x = c("CHROM", "POS", "ALT", "REF", "logp", "pval", "mincov", "maxcov", "sum_cov", "sum_mac", "maf"), by.y= c("CHROM", "POS", "ALT", "REF", "logp", "pval", "mincov", "maxcov", "sum_cov", "sum_mac", "maf"))


col_order<- c("CHROM", "POS", "TYPE", "REF", "ALT", "pval", "logp", "ANN", "mincov", "maxcov", "sum_mac", "sum_cov", "maf")
new_order2 <- c(col_order, setdiff(names(merged_sig_table), col_order))
merged_sig_table<- merged_sig_table[ , new_order2 ]
merged_sig_table<-merged_sig_table[, 1:61]
View(missing_rows)

## Now add in missing indel rows
dim(missing_rows)
missing_rows<-missing_rows[,1:59]
missing_rows$ANN<-NA
missing_rows$TYPE<-"INDEL"

new_order3 <- c(col_order, setdiff(names(missing_rows), col_order))
missing_rows2<- missing_rows[, new_order3]


# check cols are the same
identical(names(merged_sig_table), names(missing_rows2))

#bind missing rows
complete_table<-rbind(merged_sig_table, missing_rows2)
dim(complete_table)


## Sort data frame
# Your custom chromosome order
type<-c("SNP", "INDEL")
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                 "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chrmito")

# Make CHROM a factor with your desired order
complete_table$CHROM <- factor(complete_table$CHROM, levels = chrom_order)
complete_table$TYPE <- factor(complete_table$TYPE, levels = type)

# Sort the dataframe by CHROM (in your order) and then by POS (ascending)
sorted_table <- complete_table %>%
  arrange(TYPE, CHROM, POS)

dim(sorted_table)
View(sorted_table)
#------------------------------------------------------------------------------
## save table: 

write.csv(sorted_table, file="tables/SuppTable2_GenomicVariantAnn.csv", row.names=FALSE)
#------------------------------------------------------------------------------