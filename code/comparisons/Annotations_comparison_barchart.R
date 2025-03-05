########## Annotations pie chart ####################################

### load in data
dataset1<-read.csv("temp/comparisons/SNPs_annotation_counts_pie.csv")
dataset2<-read.csv("temp/comparisons/SNPs_ref_annotation_counts_pie.csv")

###

combined_data <- bind_rows(
  dataset1 %>% mutate(Dataset = "Significant Regions"),
  dataset2 %>% mutate(Dataset = "Whole Genome")
) %>%
  group_by(Dataset) %>%
  mutate(Proportion = Annotation_Count / sum(Annotation_Count))

View(combined_data)
combined_data <- combined_data %>%
  mutate(
    # Remove underscores
    Annotation = str_replace_all(Annotation, "_", " "),
    # Rename specific levels
    Annotation = recode(Annotation, "missense variant" = "nonsynonymous variant", "stop gained" = "Other", "start lost" = "Other"),
    # Relevel with the cleaned names
    Annotation = factor(Annotation, levels = c(
      "nonsynonymous variant", "synonymous variant",
      "upstream gene variant", "downstream gene variant", "Other"
    ))
  )

combined_data <- combined_data %>%
  group_by(Dataset, Annotation) %>%
  summarise(
    Annotation_Count=sum(Annotation_Count),
    Proportion= sum(Proportion),
    .groups = "drop"
  )

combined_data <- combined_data %>%
  mutate(Annotation = reorder(Annotation, Proportion)) %>%
  ungroup()

write.csv(combined_data, file="tables/annotations_table.csv")


View(combined_data)

### Barplot of proportions for comparison
barplot<-ggplot(combined_data, aes(x = Dataset, y = Proportion, fill = Annotation )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired") +  # Use the discrete viridis color scale
  labs(x = "Dataset", y = "Proportion of Annotations", fill = "Annotation") + 
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # Transparent fill for the border
    panel.grid.major = element_blank(),  # Removes major gridlines
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(
      size = 12,  # Increase font size
      angle = 0,  # Rotate the labels to avoid overlap  # Adjust horizontal justification
      vjust = 1  # Adjust vertical justification
    ),  # Increase font size of x-axis tick labels
    axis.text.y = element_text(size = 14),  # Increase font size of y-axis tick labels
    axis.title.y = element_text(size = 16),  # Increase font size of y-axis label
    axis.title.x = element_blank() 
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))  # Wrap text of x-axis labels




### I could combine more categories to make it easier to see
plot(barplot)


# Save the plot using ggsave
ggsave(filename = "figures/Figure4_annotations_barplot.pdf", plot = barplot , width = 6, height = 8)
ggsave(filename = "figures/Figure4_annotations_barplot.jpeg", plot = barplot , width = 6, height = 8)
