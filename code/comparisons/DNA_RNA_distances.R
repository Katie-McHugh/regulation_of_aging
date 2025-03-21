###########################################################################
## Create a figure to communicate distances bw variants and transcripts ##
###########################################################################

### load in data
data<- read.csv("temp_tables/RNA_comparisons_summary.csv", header=TRUE)

View(data)
## counts column gives total (e.g. the variants within 1kb are also counted 
## in the 5kb category)
## exclusive column gives just that bin (e.g. variants between 1-5 kb, but not 
## less than 1kb for the 5kb row)
#------------------------------------------------------------------------------
# Load necessary library
library(ggplot2)
library(dplyr)
#------------------------------------------------------------------------------
# Categorize upstream

### relevel bars: 
data$Position<-factor(data$Position, levels = c("genic", "upstream", "downstream"))
#data$Location<-factor(data$Location, levels = c("genic", "1kb", "5kb", "10kb", "300kb", "chrom"))
data$Location<-factor(data$Location, levels = c("chrom", "300kb", "10kb", "5kb", "1kb", "genic"))
### Barplot of proportions for comparison
barplot<-ggplot(data, aes(x = Position, y = exclusive, fill = Location)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired") +  # Use the discrete viridis color scale
  labs(x = "Position", y = "Gene Variant Count", fill = "Location") + 
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
    axis.title.y = element_text(size = 18),  # Increase font size of y-axis label
    axis.title.x = element_blank(), 
    legend.text = element_text(size = 14),  # Increase font size of legend text
    legend.title = element_text(size = 18), 
  ) +
  scale_y_continuous(breaks = seq(0, 60, by = 10))
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))  # Wrap text of x-axis labels

plot(barplot)

#ggsave(filename = "figures/Figure4_annotations_barplot_short.pdf", plot = barplot , width = 6, height = 8)
ggsave(filename = "figures/Figure5_locations_barplot.jpeg", plot = barplot , width = 6, height = 8)

