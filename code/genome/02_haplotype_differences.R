### Plot Haplotype frequencies 


## Haplotypes were estimated using code in "Haplotype_estimation_eNotebook_version" script.  Make sure to update this script (and re-read in the data below) if any changes are made to the filtering parameters.

### Haplotypes estimated 10/2/24 using file GWAS_SNPS_cov20_maf5.csv created above
founderA1=read.table ("temp/genome/founderA1_10kb_cov20_maf5.txt", header=TRUE) 
founderA2=read.table ("temp/genome/founderA2_10kb_cov20_maf5.txt", header=TRUE)
founderB3=read.table ("temp/genome/founderB3_10kb_cov20_maf5.txt", header=TRUE)
founderB4=read.table ("temp/genome/founderB4_10kb_cov20_maf5.txt", header=TRUE)

## Find the differences in haplotype frequencies between paired young and old replicates
for (i in 1:12) {
  founderA1[[paste0("diff", i)]] <- founderA1[[paste0("hap.old", sprintf("%02d", i))]] - founderA1[[paste0("hap.young", sprintf("%02d", i))]]
}

for (i in 1:12) {
  founderA2[[paste0("diff", i)]] <- founderA2[[paste0("hap.old", sprintf("%02d", i))]] - founderA2[[paste0("hap.young", sprintf("%02d", i))]]
}

for (i in 1:12) {
  founderB3[[paste0("diff", i)]] <- founderB3[[paste0("hap.old", sprintf("%02d", i))]] - founderB3[[paste0("hap.young", sprintf("%02d", i))]]
}

for (i in 1:12) {
  founderB4[[paste0("diff", i)]] <- founderB4[[paste0("hap.old", sprintf("%02d", i))]] - founderB4[[paste0("hap.young", sprintf("%02d", i))]]
}

View(founderA1)
## Find the mean differences across all pairs from the differences calculated above, and the variance
founderA1$meandiff=apply(founderA1[28:39],1,mean) # this is a first-pass simple effort to quantify consistent differences in haplotypes across replicates, per founder

founderA2$meandiff=apply(founderA2[28:39],1,mean) # this is a first-pass simple effort to quantify consistent differences in haplotypes across replicates, per founder


founderB3$meandiff=apply(founderB3[28:39],1,mean) # this is a first-pass simple effort to quantify consistent differences in haplotypes across replicates, per founder


founderB4$meandiff=apply(founderB4[28:39],1,mean) # this is a first-pass simple effort to quantify consistent differences in haplotypes across replicates, per founder

# Non-cumulative offsets #from Gaxis in Part 1B
non_cumulative_offsets <- c(0,
                            230218, 813184, 316620, 1531933, 576874, 270161, 
                            1090940, 562643, 439888, 745751, 666816, 1078177, 
                            924431, 784333, 1091291
)

# Calculate cumulative offsets
cumulative_offsets <- cumsum(non_cumulative_offsets)

# Names of the chromosomes
chromosome_names <- paste0("chr", 1:16)

# Create a named vector for cumulative offsets
named_cumulative_offsets <- setNames(cumulative_offsets, chromosome_names)
print(named_cumulative_offsets)

apply_specific_offsets <- function(df, offsets) {
  # Ensure chr is treated as a factor with the correct order
  df$chr <- factor(df$chr, levels = names(offsets))
  
  # Convert factor levels to numeric for correct indexing
  df$chr_numeric <- as.numeric(df$chr)
  
  # Assign the MB column based on cumulative offsets
  df$MB <- df$pos + offsets[levels(df$chr)][df$chr_numeric]
  
  # Drop the temporary numeric column
  df$chr_numeric <- NULL
  
  df$MB<- df$MB/1e6
  
  return(df)
}

# Apply the function to each dataset with cumulative offsets
founderA1 <- apply_specific_offsets(founderA1, named_cumulative_offsets)
founderA2 <- apply_specific_offsets(founderA2, named_cumulative_offsets)
founderB3 <- apply_specific_offsets(founderB3, named_cumulative_offsets)
founderB4 <- apply_specific_offsets(founderB4, named_cumulative_offsets)

#double check the math
numbers <- c(230218/1e6, 813184/1e6, 316620/1e6, 1531933/1e6, 576874/1e6, 270161/1e6, 
             1090940/1e6, 562643/1e6, 439888/1e6, 745751/1e6, 666816/1e6, 1078177/1e6, 
             924431/1e6, 784333/1e6, 1091291/1e6)

# Calculate the sum
total_sum <- sum(numbers)
total_sum+925347/1e6 #double check that everything aligns #12048607/1e6 is the value of the last SNP on chr16,

bottom<-0
top<-0.12

# Open the PDF device with the specified file path
pdf(file = "temp_figs/haplotype_plot.pdf", height = 10, width = 18)

par(mfrow=c(2,1))
plot(founderA1$MB,founderA1$meandiff,
     type="l",
     ylim=c(bottom,top),
     col="darkgreen",
     xlab=NA,
     xaxt= "n",
     ylab="diff in hap freq",
     main="mean diffs hap freqs old-young [1-12]")

rect(0.23,bottom,1.04,top,col="grey80",lty=0)  
rect(1.36,bottom,2.89,top,col="grey80",lty=0)  
rect(3.46,bottom,3.73,top,col="grey80",lty=0)
rect(4.82,bottom,5.39,top,col="grey80",lty=0)  
rect(5.83,bottom,6.57,top,col="grey80",lty=0)    
rect(7.24,bottom,8.32,top,col="grey80",lty=0)
rect(9.24,bottom,10.03,top,col="grey80",lty=0)
rect(11.12,bottom,12.07,top,col="grey80",lty=0) 

lines(founderA1$MB,founderA1$meandiff,col="darkgreen")
lines(founderA2$MB,founderA2$meandiff,col="red")
lines(founderB3$MB,founderB3$meandiff,col="goldenrod")
lines(founderB4$MB,founderB4$meandiff,col="blue") # these colors have been used before to refer to these specific founders

#add chromosome labels at midpts
mtext("C2",line = .5,side=1, at =0.635, cex=1.2)
mtext("C4",line = .5,side=1, at =2.125, cex=1.2)
mtext("C6",line = .5,side=1, at =3.595, cex=1.2)
mtext("C8",line = .5,side=1, at =5.105, cex=1.2)
mtext("C10",line = .5,side=1, at =6.2, cex=1.2)
mtext("C12",line = .5,side=1, at =7.78, cex=1.2)
mtext("C14",line = .5,side=1, at =9.635, cex=1.2)
mtext("C16",line = .5,side=1, at =11.595, cex=1.2)

legend("topright", 
       legend = c("DBVPG6765", "DBPVG6044","YPS128", "Y12"),    # Names of the data
       col = c("darkgreen", "red", "goldenrod", "blue"),          # Corresponding colors
       lty = c(1, 1, 1, 1),                 # Corresponding point shapes
       cex = 0.9, 
       lwd= 2,
       bg="white", 
       ncol=2)    

dev.off()
### consider branching to a version that has the variance plot? (go back to GTWAS_master_code_v2)
