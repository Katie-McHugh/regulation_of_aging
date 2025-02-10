
################################Load in Data####################################

### Load in annotations for MY data
#gene_annotation_count<-read.csv("temp/significant_annotations_count.csv") 
# includes indels
sigs<-read.csv("temp/comparisons/significant_annotations_count.csv")
ref<-read.csv("temp/comparisons/ref_annotations_count.csv")
### Load in annotations for whole genome


##### Chi-square test/ Fisher's exact test

###### Find annotations that are present in both datasets #################
common_annotations <- intersect(sigs$Annotation, ref$Annotation)
print(common_annotations)
class(common_annotations)
test_annotations<-c("downstream_gene_variant", "missense_variant", "synonymous_variant", "upstream_gene_variant")


#### Replace annotations not in common with "Other" in WG annotations
ref$Annotation <- ifelse(ref$Annotation %in% test_annotations, ref$Annotation, "Other")
sigs$Annotation <- ifelse(sigs$Annotation %in% test_annotations, sigs$Annotation, "Other")

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
#sigs <- sigs %>%
 # add_row(Annotation = "Other", Annotation_Count = 0)


### Small count values, so use Fisher's Exact Test ### too large, use simulated Fisher's exact test

dataset1<-sigs
write.csv(sigs, file= "temp/comparisons/SNPs_annotation_counts_pie.csv")
dataset2<-ref
write.csv(ref, file= "temp/comparisons/SNPs_ref_annotation_counts_pie.csv")

View(sigs)
# Create a contingency table
contingency_table <- matrix(c(
  dataset1$Annotation_Count,
  dataset2$Annotation_Count
), nrow = 2, byrow = TRUE)
View(contingency_table)

# Assign row and column names for clarity
rownames(contingency_table) <- c("Data", "Genome")
colnames(contingency_table) <- dataset1$Annotation
View(contingency_table)

fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)
##p <<0.01

chisq_result <- chisq.test(contingency_table)
print(chisq_result)

chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000)
# View results
print(fisher_result)




### Look at specific annotations

categories_of_interest <- c("upstream_gene_variant", "missense_variant", "nonsynonymous_variant")

# Subset the datasets
subset_dataset1 <- dataset1[dataset1$Annotation %in% categories_of_interest, ]
subset_dataset2 <- dataset2[dataset2$Annotation %in% categories_of_interest, ]

# Create a contingency table
contingency_table <- matrix(c(
  subset_dataset1$Annotation_Count,
  subset_dataset2$Annotation_Count
), nrow = 2, byrow = TRUE)

# Assign row and column names for clarity
rownames(contingency_table) <- c("Dataset1", "Dataset2")
colnames(contingency_table) <- subset_dataset1$Annotation # doesn't work if other value is 0

# counts of 1 or fewer grouped
combined_data <- gene_annotation_count_first %>%
  mutate(Combined_Annotations = ifelse(Annotation_Count < 2, Annotation, NA)) %>%
  mutate(Annotation = ifelse(Annotation_Count < 2, "Other", Annotation)) %>%
  group_by(Annotation) %>%
  summarise(
    Annotation_Count = sum(Annotation_Count)
  )

common_annotations <- intersect(combined_data$Annotation, gene_annotation_count_first_genome1$Annotation)
print(common_annotations)

combined_genome<-gene_annotation_count_first_genome1

# Step 2: Replace annotations not in common with "Other"
combined_genome$Annotation <- ifelse(combined_genome$Annotation %in% common_annotations, combined_genome$Annotation, "Other")

combined_genome <- combined_genome %>%
  group_by(Annotation) %>%
  summarise(Annotation_Count = sum(Annotation_Count)) %>%
  ungroup()  # Removes grouping after summarizing
View(combined_genome)


### Combined small counts, so now I can use chi-square

dataset1<-combined_data
dataset2<-combined_genome

# Create a contingency table
contingency_table <- matrix(c(
  dataset1$Annotation_Count,
  dataset2$Annotation_Count
), nrow = 2, byrow = TRUE)

# Assign row and column names for clarity
rownames(contingency_table) <- c("Data", "Genome")
colnames(contingency_table) <- dataset1$Annotation
View(contingency_table)

chisq_result <- chisq.test(contingency_table)
print(chisq_result) # does come out as significant

# Perform Fisher's Exact Test or Chi-square Test
set.seed(123)
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 1e6)

print(contingency_table)
# View results
print(fisher_result) # not working

print(gene_annotation_count_first)
### It is statistically significant, but this seems to be driven more by the categories we aren't as interested in...

