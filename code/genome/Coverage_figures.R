###############################################################################
## Supplementary Table: Coverage Boxplots across Replicates ###
###############################################################################

#------------------------------------------------------------------------------
## load in libraries
install.packages("tidyverse")
install.packages("dplyr")
library(tidyverse)


#------------------------------------------------------------------------------
## read in data

snps<-read.csv("data/clean/GWAS_SNPS_cov20_maf05.csv")
head(snps)
indels<-read.csv("data/clean/indels_cov20_maf05.csv")
head(indels)

#------------------------------------------------------------------------------
## Average coverage across all positions

snps$avg_cov<-snps$sum_cov/24

#------------------------------------------------------------------------------
## Plot

plot_coverage_all_chroms <- function(snps) {
  
  # Put chromosomes in order
  chroms <- unique(snps$CHROM)
  chrM <- "chrmito" %in% chroms
  chroms <- chroms[chroms != "chrmito"]
  chrom_nums <- as.numeric(gsub("chr", "", chroms))
  chrom_nums[is.na(chrom_nums)] <- max(chrom_nums, na.rm = TRUE) + 1 #plot chrM last
  chroms <- chroms[order(chrom_nums)]

  # Calculate offsets for each chromosome
  offsets <- c(0)
  for (i in 1:(length(chroms) - 1)) { # -1 so it stops at last chromosome
    chrom_data <- subset(snps, CHROM == chroms[i])
    offsets <- c(offsets, offsets[i] + max(chrom_data$POS) + 1) ## append to list
  }
  
  # Build a combined data frame with adjusted positions
  combined <- do.call(rbind, lapply(seq_along(chroms), function(i) {
    chrom_data <- subset(snps, CHROM == chroms[i])
    chrom_data$adj_POS <- chrom_data$POS + offsets[i] ## assign adjusted POS
    chrom_data$chrom_num <- i
    chrom_data
  }))

  # Alternate shading per chromosome
  colors <- ifelse(combined$chrom_num %% 2 == 0, "grey20", "grey45")
  
  # Base plot
  plot(
    combined$adj_POS,
    combined$avg_cov,
    pch = 20,
    cex = 0.5,
    col = colors,
    xaxt = "n",
    xlab = "Chromosome",
    ylab = "Average Coverage",
    main = "Coverage Across All Chromosomes"
  )
  
  # Add chromosome labels at the midpoint of each segment
  chrom_mids <- offsets + sapply(chroms, function(chr) {
    max(subset(snps, CHROM == chr)$POS) / 2
  })
  
  axis(1, at = chrom_mids, labels = gsub("chr", "", chroms), cex.axis = 0.7)
  
  # Add vertical dividers between chromosomes
  abline(v = offsets[-1], col = "gray80", lty = 2)

}

# Call it
plot_coverage_all_chroms(snps)

#------------------------------------------------------------------------------
## Average coverage young vs old


## calculate coverage in old reps
snps$cov_old <- rowSums(snps[, grep("^N_old", names(snps))])
snps$avg_cov_o<-snps$cov_old/12


plot_coverage_old <- function(snps) {
  
  # Put chromosomes in order
  chroms <- unique(snps$CHROM)
  chrM <- "chrmito" %in% chroms
  chroms <- chroms[chroms != "chrmito"]
  chrom_nums <- as.numeric(gsub("chr", "", chroms))
  chrom_nums[is.na(chrom_nums)] <- max(chrom_nums, na.rm = TRUE) + 1 #plot chrM last
  chroms <- chroms[order(chrom_nums)]
  
  # Calculate offsets for each chromosome
  offsets <- c(0)
  for (i in 1:(length(chroms) - 1)) { # -1 so it stops at last chromosome
    chrom_data <- subset(snps, CHROM == chroms[i])
    offsets <- c(offsets, offsets[i] + max(chrom_data$POS) + 1) ## append to list
  }
  
  # Build a combined data frame with adjusted positions
  combined <- do.call(rbind, lapply(seq_along(chroms), function(i) {
    chrom_data <- subset(snps, CHROM == chroms[i])
    chrom_data$adj_POS <- chrom_data$POS + offsets[i] ## assign adjusted POS
    chrom_data$chrom_num <- i
    chrom_data
  }))
  
  # Alternate shading per chromosome
  colors <- ifelse(combined$chrom_num %% 2 == 0, "tomato3", "salmon")
  
  # Base plot
  plot(
    combined$adj_POS,
    combined$avg_cov_o,
    pch = 20,
    cex = 0.5,
    col = colors,
    xaxt = "n",
    xlab = "Chromosome",
    ylab = "Average Coverage",
    main = "Coverage Across Old Replicates"
  )
  
  # Add chromosome labels at the midpoint of each segment
  chrom_mids <- offsets + sapply(chroms, function(chr) {
    max(subset(snps, CHROM == chr)$POS) / 2
  })
  
  axis(1, at = chrom_mids, labels = gsub("chr", "", chroms), cex.axis = 0.7)
  
  # Add vertical dividers between chromosomes
  abline(v = offsets[-1], col = "gray80", lty = 2)
  
}

plot_coverage_old(snps)
#------------------------------------------------------------------------------
## Average coverage young reps


## calculate coverage in youngreps
snps$cov_young <- rowSums(snps[, grep("^N_young", names(snps))])
snps$avg_cov_y<-snps$cov_young/12


plot_coverage_young <- function(snps) {
  
  # Put chromosomes in order
  chroms <- unique(snps$CHROM)
  chrM <- "chrmito" %in% chroms
  chroms <- chroms[chroms != "chrmito"]
  chrom_nums <- as.numeric(gsub("chr", "", chroms))
  chrom_nums[is.na(chrom_nums)] <- max(chrom_nums, na.rm = TRUE) + 1 #plot chrM last
  chroms <- chroms[order(chrom_nums)]
  
  # Calculate offsets for each chromosome
  offsets <- c(0)
  for (i in 1:(length(chroms) - 1)) { # -1 so it stops at last chromosome
    chrom_data <- subset(snps, CHROM == chroms[i])
    offsets <- c(offsets, offsets[i] + max(chrom_data$POS) + 1) ## append to list
  }
  
  # Build a combined data frame with adjusted positions
  combined <- do.call(rbind, lapply(seq_along(chroms), function(i) {
    chrom_data <- subset(snps, CHROM == chroms[i])
    chrom_data$adj_POS <- chrom_data$POS + offsets[i] ## assign adjusted POS
    chrom_data$chrom_num <- i
    chrom_data
  }))
  
  # Alternate shading per chromosome
  colors <- ifelse(combined$chrom_num %% 2 == 0, "lightblue", "steelblue4")
  
  # Base plot
  plot(
    combined$adj_POS,
    combined$avg_cov_y,
    pch = 20,
    cex = 0.5,
    col = colors,
    xaxt = "n",
    xlab = "Chromosome",
    ylab = "Average Coverage",
    main = "Coverage Across Young Replicates"
  )
  
  # Add chromosome labels at the midpoint of each segment
  chrom_mids <- offsets + sapply(chroms, function(chr) {
    max(subset(snps, CHROM == chr)$POS) / 2
  })
  
  axis(1, at = chrom_mids, labels = gsub("chr", "", chroms), cex.axis = 0.7)
  
  # Add vertical dividers between chromosomes
  abline(v = offsets[-1], col = "gray80", lty = 2)
  
}

plot_coverage_young(snps)

View(snps)

#------------------------------------------------------------------------------
## Multi-panel figure


plot_coverage_all <- function(snps) {
  
  # Set up 3 panel figure
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  
  # Sort chromosomes once, reuse for all panels
  chroms <- unique(snps$CHROM)
  chroms <- chroms[chroms != "chrmito"]
  chrom_nums <- as.numeric(gsub("chr", "", chroms))
  chrom_nums[is.na(chrom_nums)] <- max(chrom_nums, na.rm = TRUE) + 1
  chroms <- chroms[order(chrom_nums)]
  
  # Calculate offsets once, reuse for all panels
  offsets <- c(0)
  for (i in 1:(length(chroms) - 1)) {
    chrom_data <- subset(snps, CHROM == chroms[i])
    offsets <- c(offsets, offsets[i] + max(chrom_data$POS) + 1)
  }
  
  # Build combined data frame once, reuse for all panels
  combined <- do.call(rbind, lapply(seq_along(chroms), function(i) {
    chrom_data <- subset(snps, CHROM == chroms[i])
    chrom_data$adj_POS <- chrom_data$POS + offsets[i]
    chrom_data$chrom_num <- i
    chrom_data
  }))
  
  # Chromosome midpoints and labels, reused for all panels
  chrom_mids <- offsets + sapply(chroms, function(chr) {
    max(subset(snps, CHROM == chr)$POS) / 2
  })
  chrom_labels <- gsub("chr", "", chroms)
  
  # Calculate shared y axis range across all three columns
  y_max <- max(combined[, c("avg_cov", "avg_cov_o", "avg_cov_y")], na.rm = TRUE)
  y_min <- min(combined[, c("avg_cov", "avg_cov_o", "avg_cov_y")], na.rm = TRUE)
  y_range <- c(y_min, y_max)
  
  # Helper function to draw each panel
  draw_panel <- function(y_col, title, col_even, col_odd) {
    colors <- ifelse(combined$chrom_num %% 2 == 0, col_even, col_odd)
    plot(
      combined$adj_POS,
      combined[[y_col]],
      pch = 20, cex = 0.5, col = colors,
      xaxt = "n",
      ylim = y_range,
      xlab = "Chromosome",
      ylab = "Average Coverage",
      main = title
    )
    axis(1, at = chrom_mids, labels = chrom_labels, cex.axis = 0.7)
    abline(v = offsets[-1], col = "gray80", lty = 2)
  }
  
  
  # Draw the three panels
  draw_panel("avg_cov",   "Coverage Across All Chromosomes",   "grey20",   "grey45")
  draw_panel("avg_cov_o", "Coverage Across Old Replicates",    "tomato3",  "salmon")
  draw_panel("avg_cov_y", "Coverage Across Young Replicates",  "lightblue","steelblue4")
  
  # Reset plot layout
  par(mfrow = c(1, 1))
}

plot_coverage_all(snps)
