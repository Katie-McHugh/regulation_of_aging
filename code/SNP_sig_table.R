### Supplementary Table 1
### this is ONLY SNPs, not indels

## read in SNPeff annotations
ann<-read.table("data/annotated_snps.txt", header=TRUE)
snps3<-read.csv("temp/WG_CMHtest_results.csv")

## reformat ann to match working format
roman_to_chr <- c(
  "I" = "chr1", "II" = "chr2", "III" = "chr3", "IV" = "chr4", "V" = "chr5",
  "VI" = "chr6", "VII" = "chr7", "VIII" = "chr8", "IX" = "chr9", "X" = "chr10",
  "XI" = "chr11", "XII" = "chr12", "XIII" = "chr13", "XIV" = "chr14", "XV" = "chr15",
  "XVI" = "chr16", "Mito"= "chrmito"
)

# Replace Roman numerals in the CHROM column with chromosome notation
ann <- ann %>%
  mutate(CHROM = roman_to_chr[CHROM])

thresh<-(0.05/nrow(snps3))
thresh_log=-log10(thresh)

above_threshold_05 <- snps3[snps3$pval < thresh, ]
## thresholds calculated above #simplify data sets

View(snps3)

p_05<-above_threshold_05[,c("CHROM", "POS", "pval")]

p_05_ann<-merge(p_05, ann)
nrow(p_05_ann)

# Save the file to the new directory
write.csv(p_05_ann, file = "temp/sig_SNPs_p05_ann.csv", row.names = FALSE)

## refer back to these for later when looking at annotations/RNA vs DNA data
#View(p_05_ann)

