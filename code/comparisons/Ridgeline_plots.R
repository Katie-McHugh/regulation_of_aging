##########################################################
## Ridgeline Plots
##########################################################

## Based on Matt's code
## Response to reviewer comments
#---------------------------------------------------------

## Install packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.23") # Use 3.21 with R 4.5.1 or newer

install.packages("BiocManager")
BiocManager::install("rtracklayer")
library(rtracklayer)

library(remotes) #To install packages from github using 'install_github' function as below
devtools::install_github("bmansfeld/QTLseqr")
#
##Use to install bioconductor packages
BiocManager::install("QTLseqr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
BiocManager::install("MexBrewer")

# Load required packages
# From bioconductor
library(VariantAnnotation) # To parse vcf from snpeff
library(GenomicRanges) # Use to map SNP positions in results to gff information with gene ID
library(rtracklayer) # For querying with genomicranges package

# From basic CRAN
library(MexBrewer) # Colores bonito!
library(scales) # View colors
library(dplyr) # Data wrangling
library(tidyverse) # Data wrangling
library(ggpubr) # Also loads basic ggplot2
library(ggridges) # For ridgeline plots
library(cowplot) # I like the plotting theme from this package. Its clean and easy to modify.
library(patchwork) #Combine plots together 
install.packages("cowplot") 
library("cowplot")
# Set theme globally
theme_set(theme_cowplot())

#---------------------------------------------------------
########## Make summary of results from SnpEff output ###############

# Read the SnpEff-annotated VCF
# This object stores:
#  Genomic coordinates (rowRanges)
# Reference/alternate alleles (ref(), alt())
# INFO fields (including ANN) in a structured way

SNPeff_ann_long<-read.table("data/raw/annotated_snps.txt", head= TRUE) ## annotation table

# Define expected ANN fields (SnpEff spec)
ann_fields <- c(
  "Allele",
  "Effect",
  "Impact",
  "Gene_Name",
  "Gene_ID",
  "Feature_Type",
  "Feature_ID",
  "Transcript_BioType",
  "Rank",
  "HGVS_c",
  "HGVS_p",
  "cDNA_pos",
  "CDS_pos",
  "Protein_pos",
  "Distance",
  "Errors"
)


# Splits each annotation string on |
ann_matrix <- do.call(
  rbind,
  strsplit(SNPeff_ann_long$ANN, "\\|", fixed = FALSE)
)

ann_matrix <- ann_matrix[, seq_len(16), drop = FALSE]
colnames(ann_matrix) <- ann_fields

SNPeff_ann_long <- cbind(
  SNPeff_ann_long[, c("CHROM", "POS", "REF", "ALT")], # Applies meaningful column names
  as.data.frame(ann_matrix, stringsAsFactors = FALSE)
)

# Validation checks
# multiple genes per SNP? I.e, we don't leave out any cases of multiple genes affected by a single snp
any(duplicated(SNPeff_ann_long[, c("CHROM", "POS")]))

# gene IDs correct format?
head(unique(SNPeff_ann_long$Gene_ID))
#---------------------------------------------------------
all.snps.comb <- read.table("data/raw/filtered_snps.txt", header = T)

SNPeff_stats <- SNPeff_ann_long %>%
  left_join(
    all.snps.comb %>%
      dplyr::select(CHROM, POS),
    by = c("CHROM", "POS")
  )

#---------------------------------------------------------

# Katie and Molly: You may not want to do all of this filtering so you could skip or modify to filter another way if desired.

# Filter out SNPs we don't care about:
# So the ones downstream or upstream > 1000 bp
SNPeff_stats_filtered_1000 <- SNPeff_stats %>%
  filter(
    !(Effect == "downstream_gene_variant") &
      !(Effect == "upstream_gene_variant" & as.numeric(Distance) > 1000) &
      !(Effect == "intergenic_region")
  )


#---------------------------------------------------------

# Build gene-level summary
# Create list of effects labels corresponding to exonic snps
exonic_effects <- c(
  "missense_variant",
  "synonymous_variant",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "stop_retained_variant",
  "start_retained_variant",
  "missense_variant&splice_region_variant",
  "stop_gained&splice_region_variant",
  "stop_lost&splice_region_variant",
  "splice_region_variant&synonymous_variant",
  "splice_region_variant&stop_retained_variant",
  "3_prime_UTR_variant",
  "5_prime_UTR_variant",
  "5_prime_UTR_premature_start_codon_gain_variant", 
  "splice_region_variant&non_coding_transcript_exon_variant", 
  "initiator_codon_variant"
)

# Make data frame summary
gene_summary <- SNPeff_stats_filtered_1000 %>%
  group_by(Gene_ID) %>%
  summarise(
    n_SNPs = n(),
    n_in_exons = sum(Effect %in% exonic_effects, na.rm = TRUE),
    n_introns = sum(Effect %in% c(
      "intron_variant", "splice_region_variant", "splice_acceptor_variant&intron_variant",
      "splice_region_variant&intron_variant", "splice_donor_variant&intron_variant"
    ), na.rm = TRUE),
    n_upstream = sum(Effect == "upstream_gene_variant", na.rm = TRUE),
    n_modifer = sum(Impact == "MODIFIER", na.rm = TRUE),
    n_low_impact = sum(Impact == "LOW", na.rm = TRUE),
    n_mod_impact = sum(Impact == "MODERATE", na.rm = TRUE),
    n_high_impact = sum(Impact == "HIGH", na.rm = TRUE),
    effects = paste(unique(Effect), collapse = ";"),
    max_upstream_dist = ifelse(any(Effect == "upstream_gene_variant"),
                               max(as.numeric(Distance[Effect == "upstream_gene_variant"]), na.rm = TRUE),
                               NA_real_
    ),
    min_upstream_dist = ifelse(any(Effect == "upstream_gene_variant"),
                               min(as.numeric(Distance[Effect == "upstream_gene_variant"]), na.rm = TRUE),
                               NA_real_
    ),
    CHROM = unique(CHROM),
  ) %>%
  ungroup()

# After the fact, realized we should add n_in_gene column for modeling questions
gene_summary$n_in_gene <- gene_summary$n_in_exons + gene_summary$n_introns

head(gene_summary)
# Sanity check to make n_snps matches n_in_exons + n_introns + n_upstream
gene_summary_check <- gene_summary %>%
  mutate(
    classified_sum = n_in_exons + n_introns + n_upstream,
    diff = n_SNPs - classified_sum
  )

table(gene_summary_check$diff) # Should indicate 0 for all 721 genes
rm(gene_summary_check) # I don't want to save this


## # Katie and Molly: At this point, I move to using genomic to get gene lengths for my genes from the T. cali gff file
# You may not care about this block below where I calculate gene lengths and then the number of snps per kb of gene length
# But maybe its helpful so I left it in.
# Import the GFF
gff_gr <- import("data/raw/genomic.gff") ## GFF file 
# If they don't change one to fix it.
# I chose the gff to change
seqlevels(gff_gr) <- sub("^Chromosome_", "Chr_", seqlevels(gff_gr))
# Keep only gene features
genes_gr <- gff_gr[gff_gr$type == "gene"]
colnames(mcols(genes_gr))
# While we are working with the gff, grab the length of the significant genes 
#from our genomic ranges summary table
# Add gene length to GRanges metadata as a new column
mcols(genes_gr)$gene_length <- GenomicRanges::width(genes_gr)
# Pull the lengths out as a data frame
genes_lengths <- as.data.frame(genes_gr) %>%
  dplyr::select(ID, gene_length)
# Merge with gene_summary frame
gene_summary <- gene_summary %>%
  left_join(genes_lengths, by = c("Gene_ID" = "ID"))
# get rid of intermediate
rm(genes_lengths)
# sanity check, should be 0 NAs
sum(is.na(gene_summary$gene_length))
# Duplicate column in units of kb
gene_summary$gene_length_kb <- gene_summary$gene_length / 1000
# Make column with number of snps per kb of gene length
gene_summary$n_in_gene_per_kb <- gene_summary$n_in_gene / gene_summary$gene_length_kb
# Calculate ratio of snps upstream to snps in gene
gene_summary$up_to_in_ratio <- gene_summary$n_upstream / gene_summary$n_in_gene
gene_summary$up_to_in_ratio[is.infinite(gene_summary$up_to_in_ratio)] <- NA # Get rid of Inf cases where we divide by zero
gene_summary$up_to_in_ratio[is.nan(gene_summary$up_to_in_ratio)] <- NA # Get rid of Inf cases where both are zero
# Calculate ratio of snps upstream to snps in exons per kb
gene_summary$up_to_in_per_kb_ratio <- gene_summary$n_upstream / (gene_summary$n_in_exons / gene_summary$gene_length_kb)
gene_summary$up_to_in_per_kb_ratio[is.infinite(gene_summary$up_to_in_per_kb_ratio)] <- NA # Get rid of Inf cases where we divide by zero
gene_summary$up_to_in_per_kb_ratio[is.nan(gene_summary$up_to_in_per_kb_ratio)] <- NA # Get rid of Inf cases where both are zero

# Read in significant gene list from RNA seq study
sigs <- read.csv("data/clean/rnaseq_results_batch_sigs0.1_edited.csv", header = TRUE)
names(sigs)[1] <- "Gene_ID" # Rename first column
head(sigs)
gene_summary$In_sigs <- gene_summary$Gene_ID %in% sigs$Gene_ID
View(gene_summary)
gene_summary <- gene_summary[, c("Gene_ID", "n_SNPs", "CHROM", "max_upstream_dist", 
                                 "min_upstream_dist", "n_in_exons", "n_introns", 
                                 "n_upstream", "n_modifer","n_low_impact", 
                                 "n_mod_impact", "n_high_impact", "effects")]

##### Plotting and modeling #####
library(scales)
# Define color palette
pal_1 <- mex.brewer("Revolucion")
show_col(pal_1)
my_pal_ts <- pal_1[c(3, 8)]
show_col(my_pal_ts)

# Make summary data frame for overlaying mean's and CIs
# Make sqrt transformed variable for n_in_gene_per_kb in case I need it
gene_summary$n_in_gene_per_kb_tf <- sqrt(gene_summary$n_in_gene_per_kb)

# Define our groups
ridge_vars <- c(
  "n_upstream",
  "n_in_gene",
  "n_introns",
  "n_in_exons",
  "gene_length_kb",
  "n_in_gene_per_kb_tf",
  "n_modifer",
  "n_low_impact",
  "n_mod_impact",
  "n_high_impact"
)

View(gene_summary$In_sigs)
View(In_sigs)
# Pivot long to make summary
ridge_long <- gene_summary %>%
  dplyr::select(In_sigs, all_of(ridge_vars)) %>%
  pivot_longer(
    cols = all_of(ridge_vars),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.infinite(value) & !is.na(value))

View(ridge_long)

# Generate summary
ridge_summary <- ridge_long %>%
  group_by(metric, In_sigs) %>%
  summarise(
    n = sum(!is.na(value)),
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    se = sd_val / sqrt(n),
    ci_low = mean_val - qt(0.975, df = n - 1) * se,
    ci_high = mean_val + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

# Define function to auto call the correct elements when plotting
get_ridge_ci <- function(metric_name) {
  ridge_summary %>%
    filter(metric == metric_name)
}

# Upstream genes
# Run glm with quasipoisson distribution
upstream.mod <- glm(n_upstream ~ In_sigs + n_in_gene_per_kb, family = quasipoisson, data = gene_summary)
summary(upstream.mod)
exp(0.053168)

library(ggridges)

