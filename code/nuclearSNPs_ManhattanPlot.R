### Manhattan plot of SNP data from CMH test - Nuclear DNA

## Load data
snps4<-read.csv("temp/WG_CMHtest_results_nuclear.csv")

# Open the PDF device with the specified file path
pdf(file = "temp_figs/nuclearSNPs_ManhattanPlot_short.pdf", height = 5, width = 18)

x=snps4$MB
y=snps4$logp

top= 30
bottom=0

par(mar=c(7, 7, 4, 2)) 

plot(x,y,
     xlab="Position (Mb)",
     ylab="-log10p",
     main="Manhattan Plot",
     ylim=c(bottom, top),
     cex.main=2.5, #size of text for main title
     cex.lab=3, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE)  #we will do the axes by hand next

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,bottom,1.04,top,col="grey80",lty=0)  
rect(1.36,bottom,2.89,top,col="grey80",lty=0)  
rect(3.46,bottom,3.73,top,col="grey80",lty=0)
rect(4.82,bottom,5.39,top,col="grey80",lty=0)  
rect(5.83,bottom,6.57,top,col="grey80",lty=0)    
rect(7.24,bottom,8.32,top,col="grey80",lty=0)
rect(9.24,bottom,10.03,top,col="grey80",lty=0)
rect(11.12,bottom,12.07,top,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(x,y,pch=16)
box() 

#Now draw axes back in (you have more flexibility this way)
axis(1, at = c(0,12.07), labels=c(0,12.07),tick=FALSE, line= 2, cex.axis=1.5)
axis(2,cex.axis=1.5,las=2)

#add chromosome labels at midpts
mtext("C2",line = .5,side=1, at =0.635, cex=1.5)
mtext("C4",line = .5,side=1, at =2.125, cex=1.5)
mtext("C6",line = .5,side=1, at =3.595, cex=1.5)
mtext("C8",line = .5,side=1, at =5.105, cex=1.5)
mtext("C10",line = .5,side=1, at =6.2, cex=1.5)
mtext("C12",line = .5,side=1, at =7.78, cex=1.5)
mtext("C14",line = .5,side=1, at =9.635, cex=1.5)
mtext("C16",line = .5,side=1, at =11.595, cex=1.5)

### Bonferroni multiple testing correction (alpha=0.05)
thresh<-(0.05/nrow(snps3))
thresh_log=-log10(thresh)
thresh_log
abline(h = thresh_log, col = "red", lwd = 4) # alpha <0.05 #very similar to 0.05 threshold 

dev.off()
