###########################################################################
## Create a figure to communicate coding vs noncoding variants ##
###########################################################################

# Load necessary library
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)

browseVignettes("clusterProfiler")

help(clusterProfiler)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")
browseVignettes("topGO")
#------------------------------------------------------------------------------

# Data to represent the proportions
data <- data.frame(
  category = c("Blue", "Red"),
  value = c(463/801, 338/801)
)



# Create the plot
icon<-ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1, show.legend = FALSE) +
  scale_fill_manual(values = c("Blue" = "steelblue", "Red" = "salmon")) +
  coord_flip() +  # Flip coordinates to make it horizontal
  theme_void() +  # Remove axis and gridlines
  ggtitle("Horizontal Bar: 58% Blue, 42% Red")

ggsave(filename = "figures/icon_genic.jpeg", plot = icon , width = 6, height = 2)
