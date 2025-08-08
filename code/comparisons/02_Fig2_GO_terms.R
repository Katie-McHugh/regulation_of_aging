###########################################################################
## GO-analysis Visualization ##
###########################################################################

#------------------------------------------------------------------------------
## Load in packages
#------------------------------------------------------------------------------

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
#library(enrichplot)


# Install biomaRt package (if not already installed)
#iocManager::install("biomaRt")

# Load the biomaRt library
library(biomaRt)
library(ggplot2)
library(stringr)
library(stringi)
library(RColorBrewer)
library(dplyr)
library(viridis)

okabe_ito_palette <- c(
  "Black"        = "#000000",
  "Orange"       = "#E69F00",
  "sky_blue"     = "#56B4E9",
  "bluish_green" = "#009E73",
  "Yellow"       = "#F0E442",
  "Blue"         = "#0072B2",
  "vermilion"    = "#D55E00",
  "reddish_purple" = "#CC79A7"
)

#------------------------------------------------------------------------------
# load in GO-term data
#------------------------------------------------------------------------------
go_c<-read.table("results/GO_analyses/Genic/GO_genic_component.txt", 
                 header=TRUE, sep = "\t", quote = "")

go_p <- read.table("results/GO_analyses/Genic/GO_genic_process.txt", 
                   header = TRUE, sep = "\t", quote = "")

go_overlap_p<-read.table("results/GO_analyses/Overlap_List/Overlap_list_process.txt", 
                         header=TRUE, sep = "\t", quote = "")



## add column identifiers to combine lists
go_c$gene_list<-"DNA"
go_p$gene_list<-"DNA"
go_c$ann<-"component"
go_p$ann<-"process"
go_overlap_p$gene_list<- "both"
go_overlap_p$ann<-"process"


go_dna <- bind_rows(go_c, go_p, go_overlap_p)


## same plot, different colors

#if one isn't showing up, the name is probably wrong...this is how to fix it
go_dna <- go_dna %>%
  mutate(gene_list = recode(gene_list, "both" = "Shared"))                       #if one isn'y showing up, the name is probably wrong...this is how to fix it

go_dna$gene_list <- factor(go_dna$gene_list,
                           levels = c("DNA", "Shared"),
                           labels = c("DNA", "Shared"))  # You can also rename here

# Create the plot
plot2<-ggplot(go_dna, aes(x = reorder(TERM, -CORRECTED_PVALUE),
                          y = CORRECTED_PVALUE, 
                          size = NUM_LIST_ANNOTATIONS, color = gene_list)) +
  geom_point(alpha = 0.8) +     
  labs(x = "GO Terms", y = "Corrected P-value") +
  theme_minimal() +
  coord_flip() +  # Flip coordinates to make the plot horizontal
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(name = "Number of \nAnnotated Genes", 
                        range = c(3, 8))  +
  scale_color_manual(
    name = "Gene List",  # Set legend title
    values = c("DNA" =  "#56B4E9", "Shared" = "#E69F00")
  ) +
  facet_wrap(~ ann, scales = "free_y", dir = "v")  +                             #dir v stacks vertically instead of horizontally
  theme_light() +
  theme(strip.text = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.title.x = element_text(size = 16, margin = margin(t = 10)), 
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),  # Legend title size
        legend.text = element_text(size = 12))

plot(plot2)

ggsave(filename = "figures/Figure2_GO.jpeg", plot = plot2 ,
       width = 8, height = 7) 

#plot3 <- plot2 + scale_color_manual( name = "Gene List", values = brewer.pal(3, "Dark2"))

## Viridis color scheme
n <- 5
my_colors <- viridis(n + 1, option = "D")[-c(1, 3, 5)] # skip colors too light for background

plot3 <- plot2 +  scale_color_manual(values = my_colors)

#ggsave(filename = "figures/Figure2_GO_v2.jpeg", plot = plot3 , width = 8, height = 7)

### maybe try separating by process vs component vs etc with different gene 
## lists as the different colors

