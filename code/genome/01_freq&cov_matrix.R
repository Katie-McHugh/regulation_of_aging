#------------------------------------------------------------------------------
### Freq & Cov Matrix
#------------------------------------------------------------------------------

## read in data
indels<-read.csv("data/clean/indels_cov20_maf05.csv")
snps<-read.csv("data/clean/GWAS_SNPS_cov20_maf05.csv")


#------------------------------------------------------------------------------
cov<- snps %>%
   select(-starts_with("alt_"))
snps6 <- snps[, -c(3:4)]

indels_cov<- indels %>%
  select(-starts_with("alt_"))
indels2 <- indels_cov[, -c(3:4)]

## Data filtering: frequencies
freq <- apply(snps6[,3:50],1,function(x) (x[seq(1,48,2)]/(x[seq(2,49,2)]))) 
freq2 <- t(freq)  # this transposes the new matrix
freq3 <- cbind(snps6[,1:2],freq2)  
freq4<-freq3[, -c(1:2)]

cov_freq<-cbind(cov, freq4)
View(cov_freq1)
cov_freq1<-cov_freq[,-c(3:4, 29:33)]

#------------------------------------------------------------------------------
write.csv(cov_freq1, file= "results/genome/CMHtest_cov_freq_matrix.csv", row.names=FALSE)
#------------------------------------------------------------------------------

## Supplementary Table 3 : Coverage Values

### Average coverage for SNPs
View(cov)

cov_N<-cov[ , grep("^N_", names(cov))]
cov_means<-colMeans(cov_N, na.rm= TRUE)
cov_avg<-data.frame(SNPs=cov_means)
View(cov_avg)

indels_cov_N<-indels_cov[ , grep("^N_", names(indels_cov))]
indels_cov_means<-colMeans(indels_cov_N, na.rm= TRUE)
indels_cov_avg<-data.frame(indels=indels_cov_means)
View(indels_cov_avg)

cov_table<-cbind(cov_avg, indels_cov_avg)
View(cov_table)

write.table(cov_table, file="tables/ST3_GWAS_cov.txt", row.names= TRUE)

