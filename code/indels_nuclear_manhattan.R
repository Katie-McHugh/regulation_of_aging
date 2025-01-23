### Manhattan plot for nuclear indels

indels3<-read.csv("temp/indels_CMHtest_results.csv")

indels4=subset(indels3,CHROM!="chrmito")
Gaxis <- numeric(length=0)
chrs<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16')
for (k in 1:length(chrs)){
  data.samp<-subset(indels4,CHROM==chrs[k])
  if (chrs[k]=='chr1'){Gaxis.samp<-data.samp$POS}
  if (chrs[k]=='chr2'){Gaxis.samp<-data.samp$POS+230218}
  if (chrs[k]=='chr3'){Gaxis.samp<-data.samp$POS+230218+813184}
  if (chrs[k]=='chr4'){Gaxis.samp<-data.samp$POS+230218+813184+316620}
  if (chrs[k]=='chr5'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933}
  if (chrs[k]=='chr6'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874}
  if (chrs[k]=='chr7'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='chr8'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='chr9'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='chr10'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='chr11'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='chr12'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='chr13'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='chr14'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='chr15'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='chr16'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}

MB=Gaxis/1e6  ## MB stands for megabases, we are dividing by 1 million to put things on a megabase scale
indels4$MB=MB ## now we have a vector in our SNP table that we can use as a x-axis variable for plotting


### Manhattan plot for nuclear indels
nrow(indels3)
x2=indels4$MB
y2=indels4$logp

top= 30
bottom=0

# Open the PDF device with the specified file path
pdf(file = "temp_figs/indels_nuclear_manhattan.pdf", height = 10, width = 18)

par(mar=c(7, 7, 4, 2)) 

plot(x2,y2,
     xlab="Position (Mb)",
     ylab="-log10p",
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
points(x2,y2,pch=16)
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


# We talked about a combined threshold for SNPs/indels?
threshold_b_indels=0.05/nrow(indels3)#using indels3 applies the threshold to the ENTIRE genome, including the mitochondria
thresh_b_indels_log=-log10(threshold_b_indels)
abline(h = thresh_b_indels_log, col = "red", lwd = 4) #Bonferroni correction #p<0.05
dev.off()

