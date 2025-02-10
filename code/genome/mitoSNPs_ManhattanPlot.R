### Manhattan plot for mitochondrial SNPs

snps3<-read.csv("temp/genome/WG_CMHtest_results.csv")
snps_m<-subset(snps3, CHROM=="chrmito")
x2=snps_m$POS/1000
y2=snps_m$logp

top= 30
bottom=0

# Specify the file path in the output_figs directory
file_path <- file.path( "temp_figs/mitoSNPs_ManhattanPlot.pdf")

# Open the PDF device with the specified file path
pdf(file = file_path, height = 10, width = 18)

par(mar=c(7, 7, 4, 2)) 

plot(x2,y2,
     xlab="Position (kb)",
     ylab="-log10p",
     main=NA,
     ylim=c(bottom, top), #size of text for main title
     cex.axis= 1.5, 
     cey.axis=1.5,
     cex.lab=2.5, 
     cey.lab=2.5)  #we will do the axes by hand next

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(x2,y2,pch=16)
box() 

thresh<-(0.05/nrow(snps3))
thresh_log=-log10(thresh)
abline(h = thresh_log , col = "red", lwd = 4)

dev.off()