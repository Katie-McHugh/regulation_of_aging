########## Figure 5 ####################################

#------------------------------------------------------------------------------
## Load libraries
library(dplyr)
library(stringr)
library(ggplot2)

# okabe_ito_palette <- c("#D55E00",
#                        "#F0E442",
#                        "#56B4E9",
#                        "#0072B2", 
#                        "#E69F00")
#                        #"#009E73",
#                        #"#CC79A7"

okabe_ito_palette <- c("darkslategray2",
                       "darkslategray3",
                       "darkslategray4",
                       "darkslategray",
                       "black") ## change to greyscale
#"#009E73",
#"#CC79A7"

#------------------------------------------------------------------------------

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

#write.csv(combined_data, file="tables/annotations_table.csv")


### Barplot of proportions for comparison
barplot<-ggplot(combined_data, aes(x = Dataset, y = Proportion, fill = Annotation )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  # scale_fill_manual(values = colorRampPalette(c("black", "grey20", "grey40", "grey60","grey80"))(length(unique(combined_data$Annotation))),
  #                   labels = function(x) str_wrap(x, width = 12))+ 
  # scale_fill_viridis_d(
  #   option = "D",  # Options A–F; "D" is vibrant and good for categories
  #   labels = function(x) str_wrap(x, width = 12)
  # )+
  scale_fill_manual(
    values = okabe_ito_palette[1:length(unique(combined_data$Annotation))],
    labels = function(x) str_wrap(x, width = 12)
  )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +  # Wrap text of x-axis labels
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # Transparent fill for the border
    panel.grid.major = element_blank(),  # Removes major gridlines
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(
      size = 16,  # Increase font size
      angle = 0,  # Rotate the labels to avoid overlap  # Adjust horizontal justification
      vjust = 1  # Adjust vertical justification
    ),  # Increase font size of x-axis tick labels
    axis.text.y = element_text(size = 16),  # Increase font size of y-axis tick labels
    axis.title.y = element_text(size = 16),  # Increase font size of y-axis label
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 16)  # Adjust legend text size if needed
  ) +
  guides(fill = guide_legend(nrow = 5, byrow = TRUE))
?legend.position
### I could combine more categories to make it easier to see
plot(barplot)


# Save the plot using ggsave
ggsave(filename = "figures/G3_submission/Figure4_annotations_barplot_teal.tif", plot = barplot , width = 6, height = 5, dpi = 600) 
