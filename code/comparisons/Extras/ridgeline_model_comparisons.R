##########################################################
## Ridgeline Plots - Model comparisons
##########################################################

## Negative binomial vs zero inflated negative binomial
#---------------------------------------------------------
# Load required packages
# From bioconductor

#---------------------------------------------------------
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

## read in gene summary (created in Ridgeline_Fig4_panels.R)

gene_summary<-read.csv("temp/comparisons/ridgeline_gene_summary.csv", header = TRUE)

## UPSTREAM MODEL
#---------------------------------------------------------
## Negative Binomial GLM

upstream.mod <- glm.nb(n_upstream ~ In_sigs + n_in_gene_per_kb, data = gene_summary)
summary(upstream.mod)
exp(coef(upstream.mod)["In_sigsTRUE"])
exp(confint(upstream.mod)) ## large 

upstream.mod_0 <- zeroinfl(n_upstream ~ In_sigs + n_in_gene_per_kb | 1,
                         dist = "negbin", data = gene_summary)
summary(upstream.mod_0)


# Make plot
ridge_upstream <- ggplot(gene_summary, aes(x = n_upstream, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  # geom_errorbar(
  #   data = get_ridge_ci("n_upstream"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
  #   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  # ) +
  # geom_point(
  #   data = get_ridge_ci("n_upstream"), aes(x = mean_val, y = In_sigs),
  #   inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  # ) +
  labs(x = "Number of SNPs <1000 bp upstream", y = "") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_upstream))) +
  # annotate("text",
  #          x = 28, y = 2.5, size = 5,
  #          label = "OR== xx ~ ';' ~ italic(p)== xx",
  #          parse = TRUE
  # ) +
  theme(legend.position = "none")
ridge_upstream ## only six observations that aren't 0...
## not sure there's enough data here for this to be useful.


# synonymous.mod <- glm(n_synonymous ~ In_sigs + n_in_gene_per_kb, family = quasipoisson, data = gene_summary)
# #(Dispersion parameter for quasipoisson family taken to be 5.661218)

synonymous.mod <- glm.nb(n_synonymous ~ In_sigs + n_in_gene_per_kb, data = gene_summary)
summary(synonymous.mod)
exp(coef(synonymous.mod)["In_sigsTRUE"])
#p=0.004
synonymous.mod_0 <- zeroinfl(n_synonymous ~ In_sigs + n_in_gene_per_kb | 1,
                           dist = "negbin", data = gene_summary)
summary(synonymous.mod_0)
## doesn't change outcome, p = 0.003



summary(synonymous.mod)
exp(coef(synonymous.mod)["In_sigsTRUE"])

# Make plot
ridge_synonymous <- ggplot(gene_summary, aes(x = n_synonymous, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  # geom_errorbar(
  #   data = get_ridge_ci("n_synonymous"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
  #   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  # ) +
  # geom_point(
  #   data = get_ridge_ci("n_synonymous"), aes(x = mean_val, y = In_sigs),
  #   inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  # ) +
  labs(x = "Number of synonymous SNPs", y = "") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_synonymous))) +
  # annotate("text",
  #          x = 28, y = 2.5, size = 5,
  #          label = "OR== xx ~ ';' ~ italic(p)== xx",
  #          parse = TRUE
  # ) +
  theme(legend.position = "none")
ridge_synonymous



#nonsynonymous.mod <- glm(n_nonsynonymous ~ In_sigs + n_in_gene_per_kb, family = quasipoisson, data = gene_summary)
# (Dispersion parameter for quasipoisson family taken to be 4.043136)
nonsynonymous.mod <- glm.nb(n_nonsynonymous ~ In_sigs + n_in_gene_per_kb, data = gene_summary)
summary(nonsynonymous.mod) #Theta 1.4471

nonsynonymous.mod_0 <- zeroinfl(n_nonsynonymous ~ In_sigs + n_in_gene_per_kb | 1,
                             dist = "negbin", data = gene_summary)
summary(nonsynonymous.mod_0) #p = 0.06

# FROM: https://stats.oarc.ucla.edu/r/dae/negative-binomial-regression/
m3 <- glm(n_nonsynonymous ~ In_sigs + n_in_gene_per_kb, family = "poisson", data = gene_summary)
pchisq(2 * (logLik(nonsynonymous.mod) - logLik(m3)), df = 1, lower.tail = FALSE)

exp(coef(nonsynonymous.mod)["In_sigsTRUE"])

observed_zeros <- sum(gene_summary$n_nonsynonymous == 0, na.rm = TRUE)
predicted_zeros <- sum(dnbinom(0, mu = fitted(nonsynonymous.mod), size = nonsynonymous.mod$theta))
cat("Observed zeros:", observed_zeros, "\n")
cat("Predicted zeros:", round(predicted_zeros), "\n")

ridge_nonsynonymous <- ggplot(gene_summary, aes(x = n_nonsynonymous, y = In_sigs, fill = In_sigs)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height = 0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0),
    point_shape = "|", point_size = 4, point_color = "grey40"
  ) +
  # geom_errorbar(
  #   data = get_ridge_ci("n_nonsynonymous"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
  #   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  # ) +
  # geom_point(
  #   data = get_ridge_ci("n_nonsynonymous"), aes(x = mean_val, y = In_sigs),
  #   inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  # ) +
  labs(x = "Number of nonsynonymous SNPs", y = "") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$n_nonsynonymous))) +
  # annotate("text",
  #          x = 28, y = 2.5, size = 5,
  #          label = "OR== xx ~ ';' ~ italic(p)== xx",
  #          parse = TRUE
  # ) +
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
  # geom_errorbar(
  #   data = get_ridge_ci("gene_length_kb"), aes(y = In_sigs, xmin = ci_low, xmax = ci_high),
  #   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black"
  # ) +
  # geom_point(
  #   data = get_ridge_ci("gene_length_kb"), aes(x = mean_val, y = In_sigs),
  #   inherit.aes = FALSE, size = 3, shape = 21, fill = "red", color = "black"
  # ) +
  labs(x = "Gene Length in kb", y = "") +
  scale_fill_manual(values = my_pal_ts) +
  scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE" = "DE")) +
  coord_cartesian(xlim = c(0, max(gene_summary$gene_length_kb))) +
  # annotate("text",
  #          x = 75, y = 2.5, size = 5,
  #          label = "OR== xx ~ ';' ~ italic(p)== xx",
  #          parse = TRUE
  # ) +
  theme(legend.position = "none")
ridge_gene_length


ridges_patch2 <- (ridge_upstream / ridge_nonsynonymous / ridge_synonymous / ridge_gene_length) +
  plot_annotation(tag_levels = list(c("C", "D", "E", "F"))) &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 18, face = "bold", hjust = 0, vjust = 1))

ridges_patch2

### confidence intervals for coefficients: 

exp(confint(upstream.mod)) ## large 
exp(confint(synonymous.mod)) ## large 
exp(confint(nonsynonymous.mod)) ## large 


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

syn_de<-subset(syn, In_sigs == TRUE)
syn_no<-subset(syn, In_sigs == FALSE)

nrow(syn_de)
nrow(syn_no)
