### Overlaid SNPs and indels mitochondria

## load in data
indels3<-read.csv("temp/indels_CMHtest_results.csv")
indels_m=subset(indels3,CHROM=="chrmito")

snps3<-read.csv("temp/WG_CMHtest_results.csv")
snps_m=subset(snps3,CHROM=="chrmito")

# Open the PDF device with the specified file path

pdf(file = "temp_figs/indels&SNPs_mito_manhattan.pdf", height = 5, width = 10)

x1<-indels_m$POS/1000
y1<-indels_m$logp
x2 <- snps_m$POS/1000
y2 <- snps_m$logp

top <- 61
bottom <- 0

par(mar = c(4, 5, 3, 3) + 0.3)

plot(x2, y2,
     xlab = "Genomic Position (kb)",
     ylab = "-log10(p)",
     ylim = c(bottom, top),
     cex.lab = 2.2,
     xaxt = "n",
     yaxt= "n",
     type = "n",
     axes = FALSE
)

# Add x2, y2 points (SNPs)
points(x2, y2, col = "black", pch = 16)

# Add x1, y1 points (Indels) on top
points(x1, y1, col = "dodgerblue3", pch = 1)

legend("topright", 
       legend = c("SNPs", "Indels"),    # Names of the data
       col = c("black", "dodgerblue3"),          # Corresponding colors
       pch = c(16, 1),                 # Corresponding point shapes
       cex = 1, 
       bg="white")    

# Add axes
axis(1,cex.axis=1.8,las=1)
axis(2,cex.axis=1.8,las=2)

box()    # Add box around the plot

threshold_b_combined=0.05/(nrow(indels3)+nrow(snps3)) # this is a combined significance threshold
thresh_b_combined_log=-log10(threshold_b_combined)
abline(h = thresh_b_combined_log, col = "red", lwd = 4) #visually this is essentially the same as the SNP threshold

dev.off()