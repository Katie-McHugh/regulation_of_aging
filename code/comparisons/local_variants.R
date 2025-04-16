###############################################################################
## Table of potential local cis-regulatory variants
##############################################################################

## load local cis-reg data 
rna<-read.csv("results/variant_dist_table.csv")

## load annotations for significant snps
ann<-read.csv("temp/comparisons/sig_ALLannotations.csv")
View(ann)
#------------------------------------------------------------------------------

## filter for local variants
local<-subset(rna, dist_from_UGV<10000 | dist_from_DGV <10000)
local_data<-local[,c("Gene.ID", "Gene.Name", "Chromosome", "Start.Position", 
                      "End.Position", "variant_positions", 
                     "dist_from_UGV", "dist_from_DGV")]
View(local_data)
#split variant positions
library(tidyverse)
library(dplyr)

# Split rows on comma in `variant_positions`
local_clean <- local_data %>%
  separate_rows(variant_positions, sep = ",\\s*")

# Find variant locations

#split into UGV vs DGV
ugv<-subset(local_clean, dist_from_UGV<10000)

# replace NAs with variant positions for variants outside of genes
ugv_2 <- ugv %>%
  mutate(
    variant_positions = as.numeric(variant_positions),  # convert from character
    variant_positions = if_else(
      is.na(variant_positions),
      Start.Position - dist_from_UGV,
      variant_positions
    )
  )

#exclude RFA3 and WSC4, those already included in UGV list
dgv<-subset(local_clean, dist_from_DGV<10000 & dist_from_DGV>1)

# replace NAs with variant positions for variants outside of genes
dgv_2 <- dgv %>%
  mutate(
    variant_positions = as.numeric(variant_positions),  # convert from character
    variant_positions = if_else(
      is.na(variant_positions),
      Start.Position - dist_from_DGV,
      variant_positions
    )
  )

View(dgv_2)

## merge back together
variants<-rbind(ugv_2, dgv_2)

# next combine with annotations for each chrom and pos
# then look for other gene variants for those genes that are also within 10kb...
