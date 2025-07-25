#############################################################################
### Permutations for coding regions
#############################################################################

#--------------------------------------------------------------------------
## Load Libraries and functions
#--------------------------------------------------------------------------

roman_to_chr <- setNames(
  c(paste0("chr", 1:16), "chrmito"),
  c("I", "II", "III", "IV", "V", "VI",
    "VII", "VIII", "IX", "X", "XI", "XII", 
    "XIII", "XIV", "XV", "XVI", "mitochondrion"))

#--------------------------------------------------------------------------
## Load in data
#--------------------------------------------------------------------------

# Inputs
snps <- read.csv("data/clean/GWAS_SNPS_cov20_maf05.csv", header = TRUE)  # testable SNPs with coding annotation
ann<- read.table("data/raw/annotated_snps.txt", header= TRUE)
indels<-read.csv("data/clean/indels_cov20_maf05.csv", header=TRUE)
ann_indels<-read.delim("data/raw/annotated_indels.txt", head=TRUE)

#sigs_05_f<-read.csv( "temp/comparisons/sig_FIRSTannotation.csv") ### significant SNPs
#--------------------------------------------------------------------------
##Format data
#--------------------------------------------------------------------------

ann <- ann %>%
  mutate(CHROM= roman_to_chr[CHROM])

snps2<-merge(snps, ann)

snps3<-snps2[, c("CHROM", "POS", "ANN")]
View(snps3)

### Step 1: Extract the first annotation by splitting at the first comma
snps_f <- snps3 %>%
  mutate(First_ANN = sapply(strsplit(ANN, ","), `[`, 1))

### Step 2: Separate the First_ANN column by pipe character into multiple columns
snps_f <- snps_f %>%
  separate(First_ANN,
           into = c("Allele", "Annotation", "Annotation_Impact", "Gene_Name",
                    "Gene_ID", "Feature_Type", "Feature_ID",
                    "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p",
                    "cDNA.pos / cDNA.length", "CDS.pos / CDS.length",
                    "AA.pos / AA.length", "Distance", "INFO"),
           sep = "\\|", fill = "right")


ann_indels <- ann_indels %>%
  mutate(CHROM= roman_to_chr[CHROM])

indels2<-merge(indels, ann_indels)

indels3<-indels2[, c("CHROM", "POS", "ANN")]

### Step 1: Extract the first annotation by splitting at the first comma
indels_f <- indels3 %>%
  mutate(First_ANN = sapply(strsplit(ANN, ","), `[`, 1))

### Step 2: Separate the First_ANN column by pipe character into multiple columns
indels_f <- indels_f %>%
  separate(First_ANN,
           into = c("Allele", "Annotation", "Annotation_Impact", "Gene_Name",
                    "Gene_ID", "Feature_Type", "Feature_ID",
                    "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p",
                    "cDNA.pos / cDNA.length", "CDS.pos / CDS.length",
                    "AA.pos / AA.length", "Distance", "INFO"),
           sep = "\\|", fill = "right")

all_variants<-rbind(snps_f, indels_f)
nrow(all_variants)


#--------------------------------------------------------------------------
## Add genic column 
#--------------------------------------------------------------------------
ann_list<- unique(all_variants$Annotation)
print(ann_list)

# genic_annotations <- c("synonymous_variant", 
#                        "missense_variant",
#                        "stop_gained",
#                        "start_lost", 
#                        "frameshift_variant", 
#                        "frameshift_variant&start_lost", 
#                        "conservative_inframe_insertion", 
#                        "disruptive_inframe_deletion", 
#                        "disruptive_inframe_insertion" , 
#                        "frameshift_variant&stop_lost&splice_region_variant", 
#                        "conservative_inframe_deletion",  
#                        "stop_lost&splice_region_variant" ,
#                        "splice_region_variant&intron_variant", 
#                        "splice_region_variant&stop_retained_variant", 
#                        "stop_lost", 
#                        )

non_genic_ann<-c("upstream_gene_variant", "downstream_gene_variant")

all_variants <- all_variants %>%
  mutate(genic = !(Annotation %in% non_genic_ann))
all_variants<- all_variants %>%
  mutate(genic_binary = ifelse(genic, 1, 0))

View(all_variants)

#--------------------------------------------------------------------------
# all_snps$coding is TRUE/FALSE or 1/0

## Perform a permutation enrichment test to see if we are observing a higher than expected
### proportion of significant SNPs being in coding regions compared to the genome 
### on average

n_gwas <- 801  # Number of significant GWAS variants (SNPS and indels) that we identified
k_obs <-  473 # Number of GWAS SNPs in coding regions from our significant variants

# Permutation parameters
n_perm <- 1000 ### perform 100000 permutations
perm_counts <- numeric(n_perm)

set.seed(123)  # For reproducibility

## sample 801 SNPs without replacement 100000 times
for (i in 1:n_perm) {
  sampled_snps <- sample(all_variants$genic, n_gwas, replace = FALSE)
  perm_counts[i] <- sum(sampled_snps)
}

# Calculate empirical p-value
p_value <- sum(perm_counts >= k_obs) / n_perm
mean_coding_prop <- mean(all_variants$genic_binary)
expected_mean <- n_gwas * mean_coding_prop

p_value
mean_coding_prop
expected_mean

## we seem to be filtering out a lot of the noncoding regions...which probably 
## makes sense

# Visualize
hist<-hist(perm_counts, breaks = 30, main = "Permutation Test: Coding SNPs",
     xlab = "Number of coding SNPs in random set")
abline(v = k_obs, col = "red", lwd = 2)

plot(hist)
