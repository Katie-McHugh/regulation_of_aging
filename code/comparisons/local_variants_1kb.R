###############################################################################
## Table of potential local cis-regulatory variants
##############################################################################

## load local cis-reg data 
rna<-read.csv("results/variant_dist_table.csv")


## load annotations for significant snps
ann<-read.table("data/raw/annotated_snps.txt", header= TRUE)
ann_i<-read.table("data/raw/annotated_indels.txt", header=TRUE)

dna<-read.csv("data/clean/GWAS_SNPS_cov20_maf05.csv")



## convert roman numerals to chromosomes
roman_to_chr <- setNames(
  c(paste0("chr", 1:16), "chrmito"),
  c("I", "II", "III", "IV", "V", "VI",
    "VII", "VIII", "IX", "X", "XI", "XII", 
    "XIII", "XIV", "XV", "XVI", "mitochondrion"))

ann <-ann %>%                                 #fix chr format
  mutate(CHROM= roman_to_chr[CHROM])

ann_i <-ann_i %>%                                 #fix chr format
  mutate(CHROM= roman_to_chr[CHROM])


#------------------------------------------------------------------------------

## filter for local variants
local<-subset(rna, dist_from_UGV<1000 | dist_from_DGV <1000)
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
ugv<-subset(local_clean, dist_from_UGV<1000)

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
dgv<-subset(local_clean, dist_from_DGV<1000 & dist_from_DGV>1)

# replace NAs with variant positions for variants outside of genes
dgv_2 <- dgv %>%
  mutate(
    variant_positions = as.numeric(variant_positions),  # convert from character
    variant_positions = if_else(
      is.na(variant_positions),
      End.Position + dist_from_DGV,
      variant_positions
    )
  )


## merge back together
variants<-rbind(ugv_2, dgv_2)
variants$variant_positions<-as.numeric(variants$variant_positions)

# next combine with annotations for each chrom and pos
# then look for other gene variants for those genes that are also within 10kb...

## merge
ann_var<-merge(variants, ann, by.x = c("Chromosome", "variant_positions"), by.y = c("CHROM", "POS"))
ann_var_i<-merge(variants, ann_i, by.x = c("Chromosome", "variant_positions"), by.y = c("CHROM", "POS"))

all_var<-rbind(ann_var, ann_var_i)
all_var2<-all_var[,c("Chromosome", "variant_positions", "Gene.ID", "Gene.Name", "ANN")]

###############################################################################
# Look for ALL variants within 10kb of DEGS
###############################################################################

## Load in Data
#-------------------------------------------------------------------------------
sigs_05_f<-read.csv("temp/comparisons/sig_FIRSTannotation.csv")
rna_list<-read.csv("data/clean/rnaseq_results_batch_sigs0.1_edited.csv") 
## will eventually need to re-generate this file, but for efficiency I just took it from the previous run

################################################################################
#########lowerPOS is literally the smaller POS number...
#####   so lowerPOS is upstream and upperPOS is downstream!!!!!!!       ######
################################################################################
find_it_1kb <- function(data, chrom_value, lowerPOS, upperPOS) {
  # Ensure POS is numeric
  data$POS <- as.numeric(data$POS)
  
  # Filter for the specified chromosome
  filtered_data <- data[data$CHROM == chrom_value, ]
  
  # Find positions in the main gene region
  range <- filtered_data[filtered_data$POS >= lowerPOS & filtered_data$POS <= upperPOS, ]
  
  # Variants within 10kb upstream of lowerPOS
  smaller_pos <- filtered_data[filtered_data$POS >= (lowerPOS - 1000) & filtered_data$POS < lowerPOS, ]
  
  # Variants within 10kb upstream of upperPOS
  larger_pos <- filtered_data[filtered_data$POS > upperPOS & filtered_data$POS <= (upperPOS + 1000), ]
  
  # Combine results
  in_range <- if (nrow(range) > 0) paste(range$POS, collapse = "|") else NA
  upstream <- if (nrow(smaller_pos) > 0) paste(smaller_pos$POS, collapse = "|") else NA
  downstream <- if (nrow(larger_pos) > 0) paste(larger_pos$POS, collapse = "|") else NA
  
  return(list(
    in_range = in_range,
    upstream = upstream,
    downstream = downstream
  ))
}


find_sigs_3 <- function(data, range_table) {
  # Apply the `find_it` function to each row of `range_table`
  results <- apply(range_table, 1, function(row) {
    chrom_value <- row["CHROM"]
    lowerPOS <- as.numeric(row["lowerPOS"])
    upperPOS <- as.numeric(row["upperPOS"])
    
    # Call the `find_it` function
    find_it_1kb(data, chrom_value, lowerPOS, upperPOS)
  })
  
  # Convert the results (list of lists) into a data frame
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  
  # Combine the original range_table with the results data frame
  final_table <- cbind(range_table, results_df)
  
  # Return the updated table
  return(final_table)
}

resultd <- find_sigs_3(sigs_05_f, rna_list)
View(resultd)

#-------------------------------------------------------------------------------
## write table

#write.csv(as.data.frame(resultd), file = "temp_tables/local_variants_1kb.csv", row.names = FALSE, quote=FALSE)

#-------------------------------------------------------------------------------
## Great, now separate each variant into its own row again...

library(dplyr)
library(tidyr)

# Start with your full table: e.g., local_all
local_all <- resultd %>%
  pivot_longer(cols = c(in_range, upstream, downstream),
               names_to = "source",
               values_to = "variant_position") %>%
  separate_rows(variant_position, sep = "\\|\\s*") %>%
  filter(!is.na(variant_position) & variant_position != "")

local_all<-local_all[, c("Gene_ID", "Gene_Name", "CHROM", "variant_position", "source", "padj", "log2FoldChange")]

View(local_all)

#-------------------------------------------------------------------------------
## Now we can assign annotations...
View(local_all)

ann_var_all<-merge(local_all, ann, by.x = c("CHROM", "variant_position"), by.y = c("CHROM", "POS"))
ann_var_all_i<-merge(local_all, ann_i, by.x = c("CHROM", "variant_position"), by.y = c("CHROM", "POS"))

ann_local_var<-rbind(ann_var_all, ann_var_all_i)

ann_local_var<-ann_local_var[,c("Gene_ID", "Gene_Name", "CHROM", "variant_position", "source", "padj", "log2FoldChange", "ANN")]
View(ann_local_var)

#-------------------------------------------------------------------------------
## Format *another* table that describes each gene and the number of gene variants upstream, in_range, and downstream

variant_counts <- ann_local_var %>%
  group_by(Gene_Name, source) %>%
  summarise(n_variants = n(), .groups = "drop") %>%
  pivot_wider(names_from = source, values_from = n_variants, values_fill = 0)

View(variant_counts)

# Step 1: Extract the second field from the ANN column
ann_counts <- ann_local_var %>%
  mutate(ann_type = str_split(ANN, "\\|") %>% map_chr(~ .x[2])) %>%
  group_by(Gene_Name, ann_type) %>%
  summarise(n_ann = n(), .groups = "drop") %>%
  pivot_wider(names_from = ann_type, values_from = n_ann, values_fill = 0)

full_summary <- left_join(variant_counts, ann_counts, by = "Gene_Name")

View(full_summary)

variant_types <- full_summary %>%
  mutate(variant_count = in_range + upstream + downstream)

variant_types <- variant_types %>%
  mutate(variant_count = in_range + upstream + downstream) %>%
  select(-in_range, -upstream, -downstream)

View(variant_types)

#-------------------------------------------------------------------------------
## write table

#write.csv(as.data.frame(variant_types), file = "temp_tables/local_variant_types_1kb.csv", row.names = FALSE, quote=FALSE)
#write.csv(as.data.frame(full_summary), file = "temp_tables/local_variant_full_summary_1kb.csv", row.names = FALSE, quote=FALSE)
#write.csv(as.data.frame(ann_local_var), file = "temp_tables/local_variants_ann_1kb.csv", row.names = FALSE, quote=FALSE)

#-------------------------------------------------------------------------------
## summarize variant types

library(dplyr)

categorized_summary <- variant_types %>%
  # define each group by summing the columns that belong to it
  mutate(
    nonsynonymous   = missense_variant,
    synonymous      = synonymous_variant,
    inframe_indel   = rowSums(across(c(
      conservative_inframe_insertion,
      disruptive_inframe_insertion,
      disruptive_inframe_deletion
    )), na.rm = TRUE),
    intergenic      = rowSums(across(c(
      upstream_gene_variant,
      downstream_gene_variant
    )), na.rm = TRUE),
    frameshift      = frameshift_variant,
  ) %>%
  # (if you really want to drop the original per‐type columns:)
  select(-missense_variant, 
         -synonymous_variant,
         -conservative_inframe_insertion,
         -disruptive_inframe_insertion,
         -disruptive_inframe_deletion,
         -upstream_gene_variant,
         -downstream_gene_variant,
         -frameshift_variant)

#write.csv(as.data.frame(categorized_summary), file = "temp_tables/local_variants_table_1kb.csv", row.names = FALSE, quote=FALSE)

