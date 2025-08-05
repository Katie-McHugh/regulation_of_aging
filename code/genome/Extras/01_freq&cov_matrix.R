#------------------------------------------------------------------------------
### Freq & Cov Matrix
#------------------------------------------------------------------------------

## read in data
indels<-read.csv("data/clean/indels_cov20_maf05.csv")
snps<-read.csv("data/clean/GWAS_SNPS_cov20_maf05.csv")


#------------------------------------------------------------------------------
cov<- snps %>%
   select(-starts_with("alt_"))
View(cov)
snps6 <- snps[, -c(3:4)]
View(snps6)

## Data filtering: frequencies
freq <- apply(snps6[,3:50],1,function(x) (x[seq(1,48,2)]/(x[seq(2,49,2)]))) 
freq2 <- t(freq)  # this transposes the new matrix
freq3 <- cbind(snps5[,1:2],freq2)  
View(freq4)
freq4<-freq3[, -c(1:2)]

cov_freq<-cbind(cov, freq4)
View(cov_freq1)
cov_freq1<-cov_freq[,-c(3:4, 29:33)]

#------------------------------------------------------------------------------
write.csv(cov_freq1, file= "data/clean/CMHtest_cov_freq_matrix.csv", row.names=FALSE)
