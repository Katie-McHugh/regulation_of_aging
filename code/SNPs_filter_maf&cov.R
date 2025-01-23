### Filter and organize SNP table

### Load in SNP table from data folder
snps <- read.table("../data/filtered_snps.txt", header = TRUE)

snps$Nmiss=NULL

#summarize the average coverage we achieved across all the samples (so, every other column starting column 6)
coverage=seq(6,52,2) 
#summary(snps[,coverage])

## calculate the min and max coverage per site and add it as a column
snps$mincov=apply(snps[,coverage],1,min)
snps$maxcov=apply(snps[,coverage],1,max)

## filter for coverage
snps2=subset(snps,mincov>20 & maxcov<1000)  

#filter for minor allele frequency
mac=seq(5,51,2)
snps2$sum_mac=apply(snps2[,mac],1,sum)
snps2$sum_cov=apply(snps2[,coverage],1,sum)
snps2$maf=snps2$sum_mac/snps2$sum_cov

snps3=subset(snps2,maf>0.05)
snps3=subset(snps3,maf<0.95) 

snps3a=snps3
snps3a$avg_cov=as.numeric((snps3$sum_cov)/24) 
snps3a<-as.data.frame(snps3a)
summary(snps3a$avg_cov)  
avg_cov_all<-mean(snps3a$avg_cov)
print(avg_cov_all) ### average coverage across all SNPs and replicates


## if we filter for a minor allele frequency of 5%, that means we are getting rid of any SNP for which the less common allele is at a total frequency of 5% across all populations.
## importantly, this means we need to consider alleles with a total frequency (across all pops, this is the "maf" column) that is both greater than 5% and less than 95%.
## imposing this filter has the added benefit of getting rid of any non-polymorphic sites (these would be ones with a "maf" of zero or 1).  So it is more efficient than the previous version.


### Save filtered data to temp folder

# Specify the file path in the new directory
file_path1 <- file.path("temp/GWAS_SNPS_cov20_maf5.csv")

# Save the file to the new directory
write.csv(snps3, file = file_path1, row.names = FALSE)
