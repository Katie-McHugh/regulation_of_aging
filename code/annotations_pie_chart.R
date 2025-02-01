########## Annotations pie chart ####################################

### load in data
dataset1<-read.csv("temp/SNPs_annotation_counts_pie.csv")
dataset2<-read.csv("temp/SNPs_ref_annotation_counts_pie.csv")

###

combined_data <- bind_rows(
  dataset1 %>% mutate(Dataset = "Significant Regions"),
  dataset2 %>% mutate(Dataset = "Whole Genome")
) %>%
  group_by(Dataset) %>%
  mutate(Proportion = Annotation_Count / sum(Annotation_Count))


combined_data <- combined_data %>%
  mutate(
    # Remove underscores
    Annotation = str_replace_all(Annotation, "_", " "),
    # Rename specific levels
    Annotation = recode(Annotation, "missense variant" = "nonsynonymous variant"),
    # Relevel with the cleaned names
    Annotation = factor(Annotation, levels = c(
      "nonsynonymous variant", "synonymous variant",
      "upstream gene variant", "downstream gene variant", "stop gained", "start lost",
      "Other"
    ))
  )

#View(combined_data)

# Plot the two-panel pie chart
pie_chart<-ggplot(combined_data, aes(x = "", y = Proportion, fill = Annotation)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  facet_wrap(~ Dataset, strip.position = "bottom") +
  theme_void() +
  theme(legend.position = "right", 
        strip.placement = "outside",
        strip.text = element_text(size = 14, margin = margin(t = 10, b = 10))
  ) +
  scale_fill_viridis_d(option = "inferno")


pie_chart

# Save the plot using ggsave
ggsave(filename = "temp_figs/pie_chart_SNPs.pdf", plot = pie_chart, width = 8, height = 4)
