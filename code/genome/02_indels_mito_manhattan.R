### Manhattan mitochondrial indels

indels3<-read.csv("temp/genome/indels_CMHtest_results.csv")
indels_m<-subset(indels3, CHROM=="chrmito")
View(indels_m)

x2=indels_m$POS
y2=indels_m$logp

top= 30
bottom=0

# Open the PDF device with the specified file path
pdf(file = "temp_figs/indels_mito_manhattan.pdf", height = 10, width = 18)

par(mar=c(7, 7, 4, 2)) 

plot(x2,y2,
     xlab="Position (Mb)",
     ylab="-log10p",
     ylim=c(bottom, top),
     cex.main=2.5, #size of text for main title
     cex.lab=3, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE)  #we will do the axes by hand next

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(x2,y2,pch=16)
box() 
axis(1, tick=FALSE, cex.axis=1.7)
axis(2, cex.axis=1.7,las=2,)
threshold_b_indels=0.05/nrow(indels3) #using indels3 applies the threshold to the ENTIRE genome, including the mitochondria
thresh_b_indels_log=-log10(threshold_b_indels)
abline(h = thresh_b_indels_log, col = "red", lwd = 4) #Bonferroni correction #p<0.05

dev.off()