### Fisher's Test for Annotations
#------------------------------------------------------------------------------

### Load in annotations for MY data # SNPS ONLY (no indels)

## table of annotation types for significant SNPs
sigs<-read.csv("temp/comparisons/significant_annotations_count.csv")

## table of annotation types for all annotated SNPs in SNPeff
ref<-read.csv("temp/comparisons/ref_annotations_count.csv")


#------------------------------------------------------------------------------
### Reorganize datasets

###### Find annotations that are present in both datasets
common_annotations <- intersect(sigs$Annotation, ref$Annotation)
print(common_annotations)

#### Replace annotations not in common with "Other" in WG annotations
test_annotations<-c("downstream_gene_variant", "missense_variant", 
                    "synonymous_variant", "upstream_gene_variant"#, 
                    #"start_lost", "stop_gained"
                    )
head(sigs)
ref$Annotation <- ifelse(ref$Annotation %in% test_annotations, ref$Annotation, "Other")
sigs$Annotation <- ifelse(sigs$Annotation %in% test_annotations, sigs$Annotation, "Other")

head(ref)
#### reconsolidate data frame
ref<- ref %>%
  group_by(Annotation) %>%
  summarise(Annotation_Count = sum(Annotation_Count)) %>%
  ungroup()  # Removes grouping after summarizing

sigs<- sigs %>%
  group_by(Annotation) %>%
  summarise(Annotation_Count = sum(Annotation_Count)) %>%
  ungroup()  # Removes grouping after summarizing


### Add "other" column to sig data # only if it doesn't already exist
# sigs <- sigs %>%
#  add_row(Annotation = "Other", Annotation_Count = 0)

### Small count values, so use Fisher's Exact Test ### too large, use simulated Fisher's exact test

dataset1<-sigs
write.csv(sigs, file= "temp/comparisons/SNPs_annotation_counts_pie.csv")
dataset2<-ref
write.csv(ref, file= "temp/comparisons/SNPs_ref_annotation_counts_pie.csv")


#------------------------------------------------------------------------------
## Create a contingency table

# Create named vectors where names are the annotation categories
sig_data <- setNames(dataset1$Annotation_Count, dataset1$Annotation)
ref_data <- setNames(dataset2$Annotation_Count, dataset2$Annotation)

# Get all unique categories from both datasets
all_categories <- union(names(sig_data), names(ref_data))

# Ensure both tables contain all categories (fill missing ones with 0)
sig_data <- sig_data[all_categories]
ref_data <- ref_data[all_categories]

# Construct the contingency table with matching row names
contingency_table <- rbind(sig_data, ref_data)

# View the result
View(contingency_table)
#------------------------------------------------------------------------------

### Perform tests

fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)
##p <<0.01

# chisq_result <- chisq.test(contingency_table)
# print(chisq_result)

chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000)
# View results

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# ## What if smallest counts are grouped
# ref2<-ref
# sigs2<-sigs
# 
# #### Replace annotations not in common with "Other" in WG annotations
# test_annotations2<-c("downstream_gene_variant", "missense_variant", 
#                     "synonymous_variant", "upstream_gene_variant")
# 
# 
# ref2$Annotation <- ifelse(ref2$Annotation %in% test_annotations2, ref2$Annotation, "Other")
# sigs2$Annotation <- ifelse(sigs2$Annotation %in% test_annotations2, sigs2$Annotation, "Other")
# View(ref2)
# 
# #### reconsolidate data frame
# ref2<- ref2 %>%
#   group_by(Annotation) %>%
#   summarise(Annotation_Count = sum(Annotation_Count)) %>%
#   ungroup()  # Removes grouping after summarizing
# 
# sigs2<- sigs2 %>%
#   group_by(Annotation) %>%
#   summarise(Annotation_Count = sum(Annotation_Count)) %>%
#   ungroup()  # Removes grouping after summarizing
# 
# 
# ### Add "other" column to sig data # only if it doesn't already exist
# # sigs2 <- sigs2 %>%
# #   add_row(Annotation = "Other", Annotation_Count = 0)
# 
# # Create named vectors where names are the annotation categories
# sig_data2 <- setNames(dataset1$Annotation_Count, dataset1$Annotation)
# ref_data2 <- setNames(dataset2$Annotation_Count, dataset2$Annotation)
# 
# # Ensure both tables contain all categories (fill missing ones with 0)
# sig_data2 <- sig_data2[test_annotations3]
# ref_data2 <- ref_data2[test_annotations3]
# 
# # Get all unique categories from both datasets
# all_categories2 <- union(names(sig_data2), names(ref_data2))
# 
# # Construct the contingency table with matching row names
# contingency_table2 <- rbind(sig_data2, ref_data2)
# 
# # View the result
# View(contingency_table2)
# 
# fisher_result2 <- fisher.test(contingency_table2, simulate.p.value = TRUE)
# print(fisher_result2)
# ?fisher.test
# fisher_result2$estimate
# ##p <<0.01
# ## doesn't really change anything

# fisher.test(contingency_table2, workspace = 2e8)  # Increase workspace size
fisher.test(contingency_table, workspace = 2e8)  # Increase workspace size
## any way I look at it, p-value is small...

## Cramer's V for effect size? 
install.packages("effectsize")
library(effectsize)

head(contingency_table2)
head(contingency_table)
# ?cramers_v
cramers_v(contingency_table, adjust = TRUE, alternative = "two.sided")

### Post - hoc power analysis: 

##https://stats.stackexchange.com/questions/133441/computing-the-power-of-fishers-exact-test-in-r

N = sum(contingency_table)  # this is the total number of observations
N

probs = prop.table(contingency_table)       
# these are the probabilities of an observation
probs                           #  being in any given cell

probs.v = as.vector(probs)      
# notice that the probabilities read column-wise
probs.v

cuts = c(0, cumsum(probs.v))    

cuts
head(contingency_table)
set.seed(4941)      # this makes it exactly reproducible
B      = 5000    # number of iterations in simulation
vals   = runif(N*B) # generate random values / probabilities

cats   = cut(vals, breaks=cuts, labels=c("sig_DGV", "ns_DGV", "sig_MV", "ns_mv",
                                          "sig_SV", "NS_SV", "sig_UGV", 
                                         "NS_ugv", "sig_Other", "ns_Other"
                                         ))

lvls <- levels(cats)
cats   = matrix(cats, nrow=N, ncol=B, byrow=F)
counts <- apply(cats, 2, function(x) {
  as.vector(table(factor(x, levels = lvls)))
})

p.vals = vector(length=B) # this will store the outputs
ptm = proc.time()         # this lets me time the simulation
for(i in 1:B){
  mat       = matrix(counts[,i], nrow=2, ncol=5, byrow=F)
  p.vals[i] = fisher.test(mat, simulate.p.value=T)$p.value
}
proc.time() - ptm               # not too bad, really
#  user  system elapsed 
# 28.66    0.32   29.08 
#
mean(p.vals>=.05)               
# the estimated probability of type II errors is 0
# [1] 0
c(0, 3/B)                       
# using the rule of 3 to estimate the 95% CI

