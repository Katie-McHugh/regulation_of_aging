#------------------------------------------------------------------------------
### Freq & Cov Matrix
#------------------------------------------------------------------------------

## read in data
indels<-read.csv("data/clean/indels_cov20_maf05.csv")
snps<-read.csv("data/clean/GWAS_SNPS_cov20_maf05.csv")


#------------------------------------------------------------------------------
snps5<- snps %>%
  select(-starts_with("alt_"))
snps6 <- snps5[, -c(29:34)]

## Data filtering: frequencies
freq <- apply(snps3[,5:54],1,function(x) (x[seq(1,48,2)]/(x[seq(2,49,2)]))) 
freq2 <- t(freq)  # this transposes the new matrix
freq3 <- cbind(snps3[,1:2],freq2)  
View(freq3)
freq3<-freq3[, -c(1:2)]

cov_freq<-cbind(snps6, freq3)
cov_freq1<-cov_freq[,-c(30:31)]
cov_freq2<-cov_freq1[,-c(3:4)]

#------------------------------------------------------------------------------
write.csv(cov_freq2, file= "data/clean/CMHtest_cov_freq_matrix.csv", row.names=FALSE)
