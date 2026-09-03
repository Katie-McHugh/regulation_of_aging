### Site frequency spectra

## also calculate total segregating alleles

#------------------------------------------------------------------------------
# Load Data

snps<- read.csv("data/clean/GWAS_SNPS_cov20_maf05.csv")
head(snps)

## SFS Plot
#------------------------------------------------------------------------------


# pdf("figures/SFS.pdf", width = 12, height = 10)
# par(mfrow = c(6,4), mar = c(2,2,2,1))  # Adjust margins if needed


tiff(file = "figures/G3_submission/SuppFig5A.tiff", height = 6, width = 7,
     units = "in",
     res = 600)

par(mfrow = c(6,4), mar = c(2,2,2,1))
# Old replicates
hist(snps$alt_old_01/snps$N_old_01, main = "Old_01")
hist(snps$alt_old_02/snps$N_old_02, main = "Old_02")
hist(snps$alt_old_03/snps$N_old_03, main = "Old_03")
hist(snps$alt_old_04/snps$N_old_04, main = "Old_04")
hist(snps$alt_old_05/snps$N_old_05, main = "Old_05")
hist(snps$alt_old_06/snps$N_old_06, main = "Old_06")
hist(snps$alt_old_07/snps$N_old_07, main = "Old_07")
hist(snps$alt_old_08/snps$N_old_08, main = "Old_08")
hist(snps$alt_old_09/snps$N_old_09, main = "Old_09")
hist(snps$alt_old_10/snps$N_old_10, main = "Old_10")
hist(snps$alt_old_11/snps$N_old_11, main = "Old_11")
hist(snps$alt_old_12/snps$N_old_12, main = "Old_12")

## young replicates
hist(snps$alt_young_01/snps$N_young_01, main = "Young_01")
hist(snps$alt_young_02/snps$N_young_02, main = "Young_02")
hist(snps$alt_young_03/snps$N_young_03, main = "Young_03")
hist(snps$alt_young_04/snps$N_young_04, main = "Young_04")
hist(snps$alt_young_05/snps$N_young_05, main = "Young_05")
hist(snps$alt_young_06/snps$N_young_06, main = "Young_06")
hist(snps$alt_young_07/snps$N_young_07, main = "Young_07")
hist(snps$alt_young_08/snps$N_young_08, main = "Young_08")
hist(snps$alt_young_09/snps$N_young_09, main = "Young_09")
hist(snps$alt_young_10/snps$N_young_10, main = "Young_10")
hist(snps$alt_young_11/snps$N_young_11, main = "Young_11")
hist(snps$alt_young_12/snps$N_young_12, main = "Young_12")

dev.off()

## Segregating alleles
#------------------------------------------------------------------------------
head(snps)
snps$young_cov<-sum(c())

snps$young_cov <- rowSums(snps[, grepl("^N_young_", names(snps))])
snps$old_cov <- rowSums(snps[, grepl("^N_old_", names(snps))])

snps$young_alt <- rowSums(snps[, grepl("^alt_young_", names(snps))])
snps$old_alt <- rowSums(snps[, grepl("^alt_old_", names(snps))])

snps$young_seg<-snps$young_cov-snps$young_alt
snps$old_seg<-snps$old_cov-snps$old_alt

young_fixed<-subset(snps, young_seg == 0 | young_seg == young_cov)
nrow(young_fixed)

old_fixed<-subset(snps, old_seg == 0 | old_seg == old_cov)
nrow(old_fixed)


## Table of segregating alleles by sample: 
#------------------------------------------------------------------------------

# Total coverage minus alternative alleles at each site. 0s indicate fixed sites.
snps$old_1_ref<-snps$N_old_01-snps$alt_old_01
snps$old_2_ref<-snps$N_old_02-snps$alt_old_02
snps$old_3_ref<-snps$N_old_03-snps$alt_old_03
snps$old_4_ref<-snps$N_old_04-snps$alt_old_04
snps$old_5_ref<-snps$N_old_05-snps$alt_old_05
snps$old_6_ref<-snps$N_old_06-snps$alt_old_06
snps$old_7_ref<-snps$N_old_07-snps$alt_old_07
snps$old_8_ref<-snps$N_old_08-snps$alt_old_08
snps$old_9_ref<-snps$N_old_09-snps$alt_old_09
snps$old_10_ref<-snps$N_old_10-snps$alt_old_10
snps$old_11_ref<-snps$N_old_11-snps$alt_old_11
snps$old_12_ref<-snps$N_old_12-snps$alt_old_12

snps$young_1_ref<-snps$N_young_01-snps$alt_young_01
snps$young_2_ref<-snps$N_young_02-snps$alt_young_02
snps$young_3_ref<-snps$N_young_03-snps$alt_young_03
snps$young_4_ref<-snps$N_young_04-snps$alt_young_04
snps$young_5_ref<-snps$N_young_05-snps$alt_young_05
snps$young_6_ref<-snps$N_young_06-snps$alt_young_06
snps$young_7_ref<-snps$N_young_07-snps$alt_young_07
snps$young_8_ref<-snps$N_young_08-snps$alt_young_08
snps$young_9_ref<-snps$N_young_09-snps$alt_young_09
snps$young_10_ref<-snps$N_young_10-snps$alt_young_10
snps$young_11_ref<-snps$N_young_11-snps$alt_young_11
snps$young_12_ref<-snps$N_young_12-snps$alt_young_12

seg_sites<-subset(snps[, c(5, 7, 9, 11, 13, 15, 17, 19, 
                           21, 23, 25, 27, 29, 31, 33, 35, 
                           37, 39, 41, 43, 45, 47, 49, 51, 64:87)])
head(seg_sites)
## only include coverage columns and columns 


## count segregating sites: 
library(dplyr)

for (i in sprintf("%02d", 1:12)) {
  # old replicates
  seg_sites <- seg_sites %>%
    mutate("old_{i}_seg" := as.integer(
      .data[[paste0("alt_old_", i)]] != 0 &
        .data[[paste0("old_", as.integer(i), "_ref")]] != 0
    ))
  
  # young replicates
  seg_sites <- seg_sites %>%
    mutate("young_{i}_seg" := as.integer(
      .data[[paste0("alt_young_", i)]] != 0 &
        .data[[paste0("young_", as.integer(i), "_ref")]] != 0
    ))
}

seg_totals <- seg_sites %>%
  select(ends_with("_seg")) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(everything(), names_to = "replicate", values_to = "segregating_sites")

seg_totals <- seg_totals %>%
  mutate(replicate = gsub("_seg$", "", replicate))
head(seg_totals)

write.csv(seg_totals, file = "tables/supp_table_segregating_sites_by_rep.csv")
View(seg_totals)

seg_totals <- seg_totals %>%
  mutate(group = case_when(
    grepl("^young", replicate) ~ "young",
    grepl("^old",   replicate) ~ "old"
  ))

View(seg_totals)

library(ggplot2)

box_segregating<-ggplot(seg_totals, aes(x = group, y = segregating_sites, fill = group)) +
  geom_boxplot() +
  labs(x     = "Age",
       y     = "Segregating Sites"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave("figures/G3_submission/SuppFig5B.pdf", plot = box_segregating, 
       width = 4, height = 5,dpi=600)
       