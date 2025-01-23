### Load and filter indel data

ann_indels<-read.delim("data/annotated_indels.txt", head=TRUE)
indels<-read.delim("data/filtered_indels.txt", head=TRUE)


indels$Nmiss=NULL

## Use the same filtering parameters as the SNP data
#summarize the average coverage we achieved across all the samples (so, every other column starting column 6)
coverage=seq(6,52,2) 
summary(indels[,coverage])

## calculate the min and max coverage per site and add it as a column
indels$mincov=apply(indels[,coverage],1,min)
indels$maxcov=apply(indels[,coverage],1,max)
indels2=subset(indels,mincov>20 & maxcov<1000)  
nrow(indels2)

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

# Save the file to the new directory
write.csv(indels3, file = "temp/indels_cov20_maf05.csv", row.names = FALSE)

