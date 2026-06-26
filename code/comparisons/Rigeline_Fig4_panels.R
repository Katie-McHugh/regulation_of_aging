##########################################################
## Ridgeline Plots
##########################################################

## Add to Figure 4
## Based on Matt's code
## Response to reviewer comments
#---------------------------------------------------------

## Install packages

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.23") # Use 3.21 with R 4.5.1 or newer
# 
# install.packages("BiocManager")
# BiocManager::install("rtracklayer")
# #To install packages from github using 'install_github' function as below
# devtools::install_github("bmansfeld/QTLseqr")
# #
# ##Use to install bioconductor packages
# BiocManager::install("QTLseqr")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("VariantAnnotation")
# BiocManager::install("MexBrewer")
# install.packages("cowplot") 

# Load required packages
# From bioconductor

library(rtracklayer)
library(remotes) 
library(VariantAnnotation) # To parse vcf from snpeff
library(GenomicRanges) # Use to map SNP positions in results to gff information with gene ID
library(rtracklayer) # For querying with genomicranges package
library(MASS)

# From basic CRAN
library(MexBrewer) # Colores bonito!
library(scales) # View colors
library(dplyr) # Data wrangling
library(tidyverse) # Data wrangling
library(ggpubr) # Also loads basic ggplot2
library(ggridges) # For ridgeline plots
library(cowplot) # I like the plotting theme from this package. Its clean and easy to modify.
library(patchwork) #Combine plots together 
library(scales)
library(ggridges)
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

# Filter for only SNPs within 1kb of a gene:
# So the ones downstream or upstream > 1000 bp
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
head(SNPeff_stats_filtered_1000)
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
    n_downstream = sum(Effect == "downstream_gene_variant", na.rm = TRUE),
    n_nonsynonymous = sum(Effect == "missense_variant", na.rm = TRUE),
    n_synonymous = sum(Effect == "synonymous_variant", na.rm = TRUE),
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

# Sanity check to make n_snps matches n_in_exons + n_introns + n_upstream + n_downstream
gene_summary_check <- gene_summary %>%
  mutate(
    classified_sum = n_in_exons + n_introns + n_upstream + n_downstream,
    diff = n_SNPs - classified_sum
  )

table(gene_summary_check$diff) # Should indicate 0 for all 721 genes
rm(gene_summary_check) # I don't want to save this


## # Katie and Molly: At this point, I move to using genomic to get gene lengths for my genes from the T. cali gff file
# You may not care about this block below where I calculate gene lengths and then the number of snps per kb of gene length
# But maybe its helpful so I left it in.
# Import the GFF
gff_gr <- import("data/raw/genomic.gff") ## GFF file 


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


# get rid of intermediate
rm(genes_lengths)
# sanity check, should be 0 NAs
sum(is.na(gene_summary$gene_length))

# Duplicate column in units of kb
gene_summary$gene_length_kb <- gene_summary$gene_length / 1000
# Make column with number of snps per kb of gene length
gene_summary$n_in_gene_per_kb <- gene_summary$n_in_gene / gene_summary$gene_length_kb

per_kb_vals<-subset(gene_summary, 
                    #n_in_gene != "0" &
                    gene_length != "NA" &
                    n_in_gene != "NA"
                    )

nrow(gene_summary) - nrow(per_kb_vals) 
#556 removed values from warning is because they aren't within genes #okay

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

gene_summary$In_sigs <- gene_summary$Gene_ID %in% sigs$Gene_ID

##### Plotting and modeling #####

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
  "n_downstream",
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
 ## n is just counting number of total genes
## mean value gives the proportion that fall into that category...

# Define function to auto call the correct elements when plotting
get_ridge_ci <- function(metric_name) {
  ridge_summary %>%
    filter(metric == metric_name)
}

gene_summary %>%
  filter(In_sigs == TRUE) %>%
  distinct(Gene_ID)

write.csv(gene_summary, file="temp/comparisons/ridgeline_gene_summary.csv")

# Upstream genes
upstream.mod <- glm.nb(n_upstream ~ In_sigs + n_in_gene_per_kb, data = gene_summary)
summary(upstream.mod)
exp(coef(upstream.mod)["In_sigsTRUE"])

# Make plot
ridge_upstream <- ggplot(gene_summary, aes(x = n_upstream, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  labs(x = "Number of SNPs <1000 bp upstream", y = "") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_upstream))) +
  theme(legend.position = "none")
ridge_upstream ## only six observations that aren't 0...
## not sure there's enough data here for this to be useful.


######## Synonymous 

synonymous.mod <- glm.nb(n_upstream ~ In_sigs + n_in_gene_per_kb, data = gene_summary)
summary(synonymous.mod)
exp(coef(synonymous.mod)["In_sigsTRUE"])

# Make plot
ridge_synonymous <- ggplot(gene_summary, aes(x = n_synonymous, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  labs(x = "Number of synonymous SNPs", y = "") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, 
                           40))+
                           #max(gene_summary$n_synonymous))) +
  theme(legend.position = "none")
ridge_synonymous



##---------------------
#Nonsynonymous

nonsynonymous.mod <- glm.nb(n_nonsynonymous ~ In_sigs + n_in_gene_per_kb, data = gene_summary)
summary(nonsynonymous.mod)
exp(coef(nonsynonymous.mod)["In_sigsTRUE"])

ridge_nonsynonymous <- ggplot(gene_summary, aes(x = n_nonsynonymous, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  labs(x = "Number of nonsynonymous SNPs", y = "") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_nonsynonymous))) +
  theme(legend.position = "none")
ridge_nonsynonymous

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
  labs(x = "Gene Length in kb", y = "") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$gene_length_kb))) +
  theme(legend.position = "none")
ridge_gene_length


ridges_patch2 <- (ridge_upstream / ridge_nonsynonymous / ridge_synonymous / ridge_gene_length) +
  plot_annotation(tag_levels = list(c("C", "D", "E", "F"))) &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 18, face = "bold", hjust = 0, vjust = 1))

ridges_patch2


#ggsave(filename = "figures/ridgeline_fig4.png", plot = ridges_patch, width = 7, height = 12, dpi = 300)
ggsave(filename = "figures/ridgeline_fig4_panels.png", plot = ridges_patch2, width = 7, height = 12, dpi = 300)

###########################
## Fisher's Exact test ###
###########################

## upstream gene variants: 
up<-subset(gene_summary, n_upstream > 0)
up<-up[, c("Gene_ID", "n_upstream", "In_sigs")]
head(up)

up_de<-subset(up, In_sigs == TRUE)
up_no<-subset(up, In_sigs == FALSE)

nrow(up_de)
nrow(up_no)


## nonsynonymous gene variants: 
ns<-subset(gene_summary, n_nonsynonymous > 0)
ns<-ns[, c("Gene_ID", "n_nonsynonymous", "In_sigs")]
View(ns)

ns_de<-subset(ns, In_sigs == TRUE)
ns_no<-subset(ns, In_sigs == FALSE)

nrow(ns_de)
nrow(ns_no)


## synonymous gene variants: 
syn<-subset(gene_summary, n_synonymous > 0)
syn<-syn[, c("Gene_ID", "n_synonymous", "In_sigs")]
nrow(syn)

syn_de<-as.data.frame(subset(syn, In_sigs == TRUE))
syn_no<-as.data.frame(subset(syn, In_sigs == FALSE))

class(syn_de)

nrow(syn_de)
nrow(syn_no)

?geom_histogram

syn_de_hist<-ggplot(syn_de, aes(x=n_synonymous)) + geom_histogram() + scale_x_continuous(limits = c(0, 100)) 
syn_no_hist<-ggplot(syn_no, aes(x=n_synonymous)) + geom_histogram() +  scale_x_continuous(limits = c(0, 100)) 

syn_patch <- (syn_de_hist / syn_no_hist)
syn_patch
hist(syn_no$n_synonymous)

curve_de<-ggplot(syn_de, aes(x = n_synonymous)) +
  geom_density(fill = "darkslategray4", alpha = 0.6)+
  scale_x_continuous(limits = c(0, 25))

curve_no<-ggplot(syn_no, aes(x = n_synonymous)) +
  geom_density(fill = "chocolate", alpha = 0.6) +
  scale_x_continuous(limits = c(0, 25))

syn_curve<- (curve_de / curve_no)
syn_curve

