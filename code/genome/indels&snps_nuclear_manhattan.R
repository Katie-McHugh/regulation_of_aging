### Plot of nuclear SNPs/indels overlaid

### !!! need to have MB scale for read-in data

## load in snp and indel data
indels3<-read.csv("temp/genome/indels_CMHtest_results.csv")
indels4<-read.csv("temp/genome/indels_CMHtest_results_nuclear.csv")
snps3<-read.csv("temp/genome/WG_CMHtest_results.csv")
snps4<-read.csv("temp/genome/WG_CMHtest_results_nuclear.csv")

View(snps4)

x1<-indels4$MB
y1<-indels4$logp
x2 <- snps4$MB
y2 <- snps4$logp

top <- 30
bottom <- 0

# Open the PDF device with the specified file path
pdf(file = "figures/indels&SNPs_nuclear_manhattan.pdf", height = 5, width = 10)

par(mar = c(4, 5, 3, 3) + 0.3)

plot(x2, y2,
     xlab = "Genomic Position",
     ylab = "-log10(p)",
     ylim = c(bottom, top),
     cex.lab = 2,
     xaxt = "n",
     yaxt= "n",
     type = "n",
     axes = FALSE
)

# Draw gray rectangles that delineate chromosomes
rect(0.23,bottom,1.04,top,col="grey92",lty=0)  
rect(1.36,bottom,2.89,top,col="grey92",lty=0)  
rect(3.46,bottom,3.73,top,col="grey92",lty=0)
rect(4.82,bottom,5.39,top,col="grey92",lty=0)  
rect(5.83,bottom,6.57,top,col="grey92",lty=0)    
rect(7.24,bottom,8.32,top,col="grey92",lty=0)
rect(9.24,bottom,10.03,top,col="grey92",lty=0)
rect(11.12,bottom,12.07,top,col="grey92",lty=0) 

# Add chromosome labels at midpts
mtext("C2",line = .5,side=1, at =0.635, cex=1.2, col="grey20")
mtext("C4",line = .5,side=1, at =2.125, cex=1.2, col="grey20")
mtext("C6",line = .5,side=1, at =3.595, cex=1.2, col="grey20")
mtext("C8",line = .5,side=1, at =5.105, cex=1.2, col="grey20")
mtext("C10",line = .5,side=1, at =6.2, cex=1.2, col="grey20")
mtext("C12",line = .5,side=1, at =7.78, cex=1.2, col="grey20")
mtext("C14",line = .5,side=1, at =9.635, cex=1.2, col="grey20")
mtext("C16",line = .5,side=1, at =11.595, cex=1.2, col="grey20")


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
axis(2,cex.axis=1.5,las=2)
box()    # Add box around the plot

### draw in significance threshold
threshold_b_combined=0.05/(nrow(indels3)+nrow(snps3)) # this is a combined significance threshold
thresh_b_combined_log=-log10(threshold_b_combined)
abline(h = thresh_b_combined_log, col = "red", lwd = 4) #visually this is essentially the same as the SNP threshold

dev.off()