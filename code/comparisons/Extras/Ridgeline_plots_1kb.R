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
#To install packages from github using 'install_github' function as below
devtools::install_github("bmansfeld/QTLseqr")
#
##Use to install bioconductor packages
BiocManager::install("QTLseqr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
BiocManager::install("MexBrewer")
install.packages("cowplot") 

# Load required packages
# From bioconductor

library(rtracklayer)

library(remotes) 
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
# So the ones downstream or upstream > 5000 bp
SNPeff_stats_filtered_1000 <- SNPeff_stats %>%
  filter(
    !(Effect == "downstream_gene_variant" & as.numeric(Distance) > 1000) &
      !(Effect == "upstream_gene_variant" & as.numeric(Distance) > 1000) &
      !(Effect == "intergenic_region" & as.numeric(Distance) > 1000)
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
View(gff_gr$ID)

head(gene_summary$Gene_ID)



# If they don't match change one to fix it.
# I chose the gff to change

###DOUBLE CHECK THESE####
chrom_map <- c(
  "NC_001133.9"  = "I",
  "NC_001134.8"  = "II",
  "NC_001135.5"  = "III",
  "NC_001136.10" = "IV",
  "NC_001137.3"  = "V",
  "NC_001138.5"  = "VI",
  "NC_001139.9"  = "VII",
  "NC_001140.6"  = "VIII",
  "NC_001141.2"  = "IX",
  "NC_001142.9"  = "X",
  "NC_001143.9"  = "XI",
  "NC_001144.5"  = "XII",
  "NC_001145.3"  = "XIII",
  "NC_001146.8"  = "XIV",
  "NC_001147.6"  = "XV",
  "NC_001148.4"  = "XVI",
  "NC_001224.1"  = "Mito"
)

seqlevels(gff_gr) <- chrom_map[seqlevels(gff_gr)]

# Keep only gene features
genes_gr <- gff_gr[gff_gr$type == "gene"]
seqlevels(gff_gr)
colnames(mcols(genes_gr))
# While we are working with the gff, grab the length of the significant genes 
#from our genomic ranges summary table
# Add gene length to GRanges metadata as a new column
mcols(genes_gr)$gene_length <- GenomicRanges::width(genes_gr)
head(genes_gr)
# Pull the lengths out as a data frame
genes_lengths <- as.data.frame(genes_gr) %>%
  dplyr::select(ID, gene_length)

genes_lengths <- genes_lengths %>%
  mutate(ID = sub("^gene[:-]", "", ID))

# Merge with gene_summary frame
gene_summary <- gene_summary %>%
  left_join(genes_lengths, by = c("Gene_ID" = "ID"))

View(gene_summary)
# get rid of intermediate
rm(genes_lengths)
# sanity check, should be 0 NAs
sum(is.na(gene_summary$gene_length))

View(gene_summary)
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
View(gene_summary)
# Read in significant gene list from RNA seq study
sigs <- read.csv("data/clean/rnaseq_results_batch_sigs0.1_edited.csv", header = TRUE)
names(sigs)[1] <- "Gene_ID" # Rename first column
head(sigs)
gene_summary$In_sigs <- gene_summary$Gene_ID %in% sigs$Gene_ID
View(gene_summary)
#gene_summary <- gene_summary[, c("Gene_ID", "In_sigs", "n_SNPs", "CHROM", "max_upstream_dist", 
#                                 "min_upstream_dist", "n_in_exons", "n_introns", 
#                                "n_upstream", "n_modifer","n_low_impact", 
#                                "n_mod_impact", "n_high_impact", "n_in_gene", "effects")]

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
colnames(gene_summary)
# Define our groups
ridge_vars <- c(
  "n_upstream",
  "n_in_gene",
  "n_introns",
  "n_in_exons",
  "gene_length_kb",
  "n_in_gene_per_kb",
  "n_modifer",
  "n_low_impact",
  "n_mod_impact",
  "n_high_impact"
)

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
View(gene_summary)

head(ridge_long)


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
View(ridge_summary) ## n is just counting number of total genes
## mean value gives the proportion that fall into that category...

# Define function to auto call the correct elements when plotting
get_ridge_ci <- function(metric_name) {
  ridge_summary %>%
    filter(metric == metric_name)
}

gene_summary %>%
  filter(In_sigs == TRUE) %>%
  distinct(Gene_ID)

View(gene_summary)

# Upstream genes
# Run glm with quasipoisson distribution
upstream.mod <- glm(n_upstream ~ In_sigs + n_in_gene_per_kb, family = quasipoisson, data = gene_summary)
summary(upstream.mod)
exp(coef(upstream.mod)["In_sigsTRUE"])

library(ggridges)


# Make plot
ridge_upstream <- ggplot(gene_summary, aes(x = n_upstream, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  geom_errorbar(
    data = get_ridge_ci("n_upstream"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
    orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  ) +
  geom_point(
    data = get_ridge_ci("n_upstream"), aes(x = mean_val, y = In_sigs),
    inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  ) +
  labs(x = "Number of SNPs 1000 bp upstream", y = "", title = "Distribution of Upstream SNP Counts") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_upstream))) +
  annotate("text",
           x = 28, y = 2.5, size = 5,
           label = "OR== xx ~ ';' ~ italic(p)== xx",
           parse = TRUE
  ) +
  theme(legend.position = "none")
ridge_upstream


# In gene
# Run glm with poisson distribution
ingene.mod <- glm(n_in_gene ~ In_sigs + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
summary(ingene.mod)
exp(coef(ingene.mod)["In_sigsTRUE"])

# Make plot
ridge_in_gene <- ggplot(gene_summary, aes(x = n_in_gene, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 5, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  geom_errorbar(
    data = get_ridge_ci("n_in_gene"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
    orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  ) +
  geom_point(
    data = get_ridge_ci("n_in_gene"), aes(x = mean_val, y = In_sigs),
    inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  ) +
  labs(x = "Number of SNPs in Genes", y = "", title = "Distribution of SNP counts in genes") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_in_gene))) +
  annotate("text",
           x = 250, y = 2.5, size = 5,
           label = "OR==1.21 ~ ';' ~ italic(p)<0.001",
           parse = TRUE
  ) +

  theme(legend.position = "none")
ridge_in_gene

# In gene per kb
ingeneperkb.mod <- lm(n_in_gene_per_kb_tf ~ In_sigs, data = gene_summary)
summary(ingeneperkb.mod)
plot(ingeneperkb.mod)
# Make plot
ridge_in_gene_per_kb <- ggplot(gene_summary, aes(x = n_in_gene_per_kb_tf, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 0.2, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  geom_errorbar(
    data = get_ridge_ci("n_in_gene_per_kb_tf"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
    orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  ) +
  geom_point(
    data = get_ridge_ci("n_in_gene_per_kb_tf"), aes(x = mean_val, y = In_sigs),
    inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  ) +
  labs(x = "Number of SNPs in genes per kb", y = "", title = "Distribution of SNP counts in genes per kb") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_in_gene_per_kb_tf))) +
  annotate("text",
           x = 4.4, y = 2.9, size = 5,
           label = "italic(beta)==0.29 ~ ';' ~ italic(p)==0.003",
           parse = TRUE
  ) +

  theme(legend.position = "none")
ridge_in_gene_per_kb

# Gene length
# Run glm with poisson distribution
length.mod <- glm(gene_length_kb ~ In_sigs, family = Gamma(link = "log"), data = gene_summary)
summary(length.mod)
exp(coef(length.mod)["In_sigsTRUE"])
# Make plot
ridge_gene_length <- ggplot(gene_summary, aes(x = gene_length_kb, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 5, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  geom_errorbar(
    data = get_ridge_ci("gene_length_kb"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
    orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  ) +
  geom_point(
    data = get_ridge_ci("gene_length_kb"), aes(x = mean_val, y = In_sigs),
    inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  ) +
  labs(x = "Gene Length in kb", y = "", title = "Distribution of Gene Lengths") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$gene_length_kb))) +
  annotate("text",
           x = 75, y = 2.5, size = 5,
           label = "OR== xx ~ ';' ~ italic(p)== xx",
           parse = TRUE
  ) +
 theme(legend.position = "none")
ridge_gene_length ## Warning message:Removed 556 rows containing non-finite outside the scale range (`stat_density_ridges()`). 
View(gene_summary[!is.finite(gene_summary$gene_length_kb), ])
# In exons
# Run glm with poisson distribution
inexon.mod <- glm(n_in_exons ~ In_sigs + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
summary(inexon.mod)
exp(coef(inexon.mod)["In_sigsTRUE"])

# Make plot
ridge_in_exon <- ggplot(gene_summary, aes(x = n_in_exons, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  geom_errorbar(
    data = get_ridge_ci("n_in_exons"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
    orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  ) +
  geom_point(
    data = get_ridge_ci("n_in_exons"), aes(x = mean_val, y = In_sigs),
    inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  ) +
  labs(x = "# SNPs in exons", y = "", title = "Distribution of SNP counts in Exons") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_in_exons))) +
  annotate("text",
           x = 70, y = 2.5, size = 5,
           label = "OR== xx ~ ';' ~ italic(p)== xx ",
           parse = TRUE
  ) +
  theme(legend.position = "none")
ridge_in_exon

# In introns
# Run glm with poisson distribution
inintron.mod <- glm(n_introns ~ In_sigs + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
summary(inintron.mod)
exp(coef(inintron.mod)["In_sigsTRUE"])

# Make plot
ridge_in_intron <- ggplot(gene_summary, aes(x = n_introns, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 4, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  geom_errorbar(
    data = get_ridge_ci("n_introns"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
    orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  ) +
  geom_point(
    data = get_ridge_ci("n_introns"), aes(x = mean_val, y = In_sigs),
    inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  ) +
  labs(x = "# SNPs in introns", y = "", title = "Distribution of SNP counts in Introns") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_introns))) +
  annotate("text",
           x = 220, y = 2.5, size = 5,
           label = "OR== xx ~ ';' ~ italic(p)< xx",
           parse = TRUE
  ) +
  theme(legend.position = "none")
ridge_in_intron


ridges_patch <- (ridge_in_gene_per_kb / ridge_in_exon / ridge_in_intron / ridge_upstream / ridge_gene_length)

ggsave(filename = "figures/ridgeline_column_plot_1kb.png", plot = ridges_patch, width = 7, height = 12, dpi = 300)


