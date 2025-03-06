
################################################################################
# Load and filter indel data
################################################################################

## SNPeff annotation file
ann_indels<-read.delim("data/raw/annotated_indels.txt", head=TRUE)

## indel table
indels<-read.delim("data/raw/filtered_indels.txt", head=TRUE)

#------------------------------------------------------------------------------

## remove Nmiss column 
indels$Nmiss=NULL                                 ## already filtered, all= 0  

## Use same filters as SNP data

coverage=seq(6,52,2)                              ## every other col, start at 6
summary(indels[,coverage])                        ## summarize coverage


indels$mincov=apply(indels[,coverage],1,min)      ## add min cov col
indels$maxcov=apply(indels[,coverage],1,max)      ## add max cov col
indels2=subset(indels,mincov>20 & maxcov<1000)    ## filter for cov
nrow(indels2) 

#------------------------------------------------------------------------------

#now we can filter for minor allele frequency (we talked about this verbally but didn't do it yet)
mac=seq(5,51,2)
indels2$sum_mac=apply(indels2[,mac],1,sum)
indels2$sum_cov=apply(indels2[,coverage],1,sum)
indels2$maf=indels2$sum_mac/indels2$sum_cov

indels3=subset(indels2,maf>0.05) # this is only 87 SNPs
indels3=subset(indels3,maf<0.95)

indels3a=indels3
indels3a$avg_cov=as.numeric((indels3$sum_cov)/24) 
indels3a<-as.data.frame(indels3a)
summary(indels3a$avg_cov)  
avg_cov_all_indels<-mean(indels3a$avg_cov)
print(avg_cov_all_indels) ### average coverage across all indels and replicates

#------------------------------------------------------------------------------

# Save the file to the new directory
write.csv(indels3, file = "data/clean/indels_cov20_maf05.csv", row.names = FALSE)

