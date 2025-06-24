###############################################################################
## Coverage and Frequency Plots
###############################################################################

# Libraries
#------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(devtools)
library(poolSeq)
#install.packages("cowplot")
library(cowplot)
#install.packages("patchwork")
library(patchwork)
#------------------------------------------------------------------------------

# Read in data
#------------------------------------------------------------------------------
cf =read.csv("data/clean/CMHtest_cov_freq_matrix.csv", header= TRUE)
data<-read.table("data/raw/annotated_indels.txt")
View(cf)
#------------------------------------------------------------------------------


# Creating coverage matrix
#------------------------------------------------------------------------------
coverage<-cf[, c(1:27)]

c_long <- reshape2::melt(coverage, id.vars = c("CHROM", "POS", "logp"))
cov_long <- c_long %>%
  mutate(age = ifelse(grepl("old", variable), "old", "young"),
         pair = as.numeric(gsub("[^0-9]", "", variable)))
#------------------------------------------------------------------------------

# Frequency matrix
#------------------------------------------------------------------------------
freq<-cf[, c(1:2,27,30:53)]


f_long <- reshape2::melt(freq, id.vars = c("CHROM", "POS", "logp"))
freq_long <- f_long %>%
  mutate(age = ifelse(grepl("old", variable), "old", "young"),
         pair = as.numeric(gsub("[^0-9]", "", variable)))

View(freq_long)
#------------------------------------------------------------------------------

# Frequency Plot
#------------------------------------------------------------------------------
#set SNP
SNP<-59253
chrom<-"chr11"

#Create Plot
freq_snp <- freq_long %>% filter(POS == SNP, CHROM == chrom) 
cov_snp <- cov_long %>% filter(POS == SNP, CHROM == chrom) 

print(class(freq_long$CHROM))
print(class(chrom))

#Create Plot
freq_plot<-ggplot(freq_snp, aes(x = factor(age, levels= c("young", "old")), 
                                y = value, group = pair))+
  geom_point(aes(color = age), size = 3)+
  geom_line()+
  labs(y = "Allele Frequency")+
  scale_x_discrete(labels = c("young" = "young", "old" = "aged")) +
  theme(strip.placement = "outside",
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))

cov_plot<-ggplot(cov_snp, aes(x = factor(age, levels = c("young", "old")), 
                              y = value, group = pair))+
  geom_point(aes(color = age), size = 3)+
  geom_line()+
  labs(y = "Coverage")+
  scale_x_discrete(labels = c("young" = "young", "old" = "aged")) +
  theme(strip.placement = "outside",
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))

merged_plot <- plot_grid(freq_plot, cov_plot, ncol = 2)


ggsave("chr11pos68563_tor2.pdf", plot = merged_plot, height = 4, width = 8)

print(merged_plot)
#------------------------------------------------------------------------------

#freq plot 2
#------------------------------------------------------------------------------
SNP<-59253
chrom<-"chr11"

freq_snp <- freq_long %>% filter(POS == SNP, CHROM == chrom) 
cov_snp <- cov_long %>% filter(POS == SNP, CHROM == chrom) 

# Create Plot
freq_plot <- ggplot(freq_snp, aes(x = factor(age, levels = c("young", "old")), 
                                  y = value, group = pair)) +
  geom_point(aes(color = age), size = 3) +
  geom_line() +
  labs(y = "Allele Frequency") +
  scale_x_discrete(labels = c("young" = "young", "old" = "aged")) +
  scale_color_manual(values = c("young" = "salmon", "old" = "steelblue")) +
  scale_y_continuous(
    limits = c(0, 1)  # Set your desired y-axis limits here
  ) +
  theme(
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 18),  # Adjust font size
    axis.title.y = element_text(size = 18),  # Adjust font size of y-axis title
    panel.background = element_rect(fill="white"),  # Remove panel background
    plot.background = element_rect(fill="white")  # Remove plot background
  )

print(freq_plot)
ggsave("chr11pos96379_tor2.pdf", plot = freq_plot, height = 3, width = 2.5)
library(svglite)

svglite::svglite("figures/tor2.svg", width = 3, height = 4)
freq_plot
dev.off()

SNP<-68563
chrom<-"chr11"

freq_snp <- freq_long %>% filter(POS == SNP, CHROM == chrom) 
cov_snp <- cov_long %>% filter(POS == SNP, CHROM == chrom) 

# Create Plot
freq_plot <- ggplot(freq_snp, aes(x = factor(age, levels = c("young", "old")), 
                                  y = value, group = pair)) +
  geom_point(aes(color = age), size = 3) +
  geom_line() +
  labs(y = "Allele Frequency") +
  scale_x_discrete(labels = c("young" = "young", "old" = "aged")) +
  scale_color_manual(values = c("young" = "salmon", "old" = "steelblue")) +
  scale_y_continuous(
    limits = c(0, 1)  # Set your desired y-axis limits here
  ) +
  theme(
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 16),  # Adjust font size
    axis.title.y = element_text(size = 16),  # Adjust font size of y-axis title
    panel.background = element_rect(fill="white"),  # Remove panel background
    plot.background = element_rect(fill="white")  # Remove plot background
  )

print(freq_plot)
ggsave("chr11pos96379_ptk1.pdf", plot = freq_plot, height = 3, width = 2.5)

cov_plot <- ggplot(cov_snp, aes(x = factor(age, levels = c("young", "old")), 
                                y = value, group = pair)) +
  geom_point(aes(color = age), size = 3) +
  geom_line() +
  labs(y = "Coverage") +
  scale_x_discrete(labels = c("young" = "young", "old" = "aged")) +
  theme(
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 16),  # Adjust font size
    axis.title.y = element_text(size = 16),  # Adjust font size of y-axis title
    panel.background = element_blank(),  # Remove panel background
    plot.background = element_blank()  # Remove plot background
  )

merged_plot <- plot_grid(freq_plot, cov_plot, ncol = 2)

ggsave("chr11pos59253_freq.pdf", plot = merged_plot, height = 4, width = 6)
#------------------------------------------------------------------------------
