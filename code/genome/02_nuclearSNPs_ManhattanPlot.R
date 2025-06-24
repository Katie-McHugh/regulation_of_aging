### Manhattan plot of SNP data from CMH test - Nuclear DNA
#-------------------------------------------------------------------------------

## Load data
snps4<-read.csv("temp/genome/WG_CMHtest_results_nuclear.csv")

#-------------------------------------------------------------------------------
 # add FWER and FDR corrected pvalues to snps4

snps4$fdr<-p.adjust(snps4$pval, method="fdr")
snps4$logFDR = (-log10(snps4$fdr))
snps4$fwer<-p.adjust(snps4$pval, method="bonferroni")
snps4$logFWER = (-log10(snps4$fwer))

#-------------------------------------------------------------------------------
# Manhattan plot of raw p-values with bonferroni corrected threshold

# Open the PDF device with the specified file path
pdf(file = "temp_figs/nuclearSNPs_ManhattanPlot.pdf", height = 5, width = 18)
?pdf()

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
bthresh<-(0.05/nrow(snps4))
bthresh_log=-log10(bthresh)
bthresh_log
abline(h = bthresh_log, col = "red", lwd = 4) # alpha <0.05 #very similar to 0.05 threshold 

dev.off()

#-------------------------------------------------------------------------------
# Manhattan plot of FWER (bonferroni) corrected p-values with a=0.05 threshold

# Open the PDF device with the specified file path
pdf(file = "temp_figs/nuclearSNPs_ManhattanPlot_FWER.pdf", height = 5, width = 18)

x=snps4$MB
y=snps4$logFWER

top= 20
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
thresh<-0.05
thresh_log=-log10(thresh)
thresh_log
abline(h = thresh_log, col = "red", lwd = 4) # alpha <0.05 #very similar to 0.05 threshold 

dev.off()

View(snps4)

#-------------------------------------------------------------------------------
# Manhattan plot of FDR (benjamini-hochberg) corrected p-values with a=0.05 threshold

# Open the PDF device with the specified file path
pdf(file = "temp_figs/nuclearSNPs_ManhattanPlot_FDR.pdf", height = 5, width = 18)

x=snps4$MB
y=snps4$logFDR

top= 20
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
thresh<-0.05
thresh_log=-log10(thresh)
thresh_log

thresh2<-0.01
thresh2_log=-log10(thresh2)
thresh2_log
abline(h = thresh_log, col = "red", lwd = 4) # alpha <0.05 #very similar to 0.05 threshold 
abline(h = thresh2_log, col = "blue", lwd = 4) # alpha <0.05 #very similar to 0.05 threshold 

dev.off()

View(snps4)

#-------------------------------------------------------------------------------
## Subset significant data

## using FDR correction on p-values
sig_snps<-subset(snps4, logFWER > thresh_log)

View(snps4)
## using FDR correction on threshold
sig_snps1<-subset(snps4, logp > bthresh_log)

nrow(sig_snps1)
nrow(sig_snps)

View(sig_snps1)
### find any differences between two methods: 

library(dplyr)

diff_1 <- anti_join(sig_snps1, sig_snps)
diff_2 <- anti_join(sig_snps, sig_snps1)

# Rows in df1 but not in df2
print(diff_1)

# Rows in df2 but not in df1
print(diff_2)

all_differences <- bind_rows(diff_1, diff_2)
print(all_differences) # two data frames are identical #yay!

# method 2
comparison <- sig_snps != sig_snps1
which(comparison, arr.ind = TRUE)

## method 3
install.packages("waldo")  # If not already installed
library(waldo)

compare(sig_snps, sig_snps1)
