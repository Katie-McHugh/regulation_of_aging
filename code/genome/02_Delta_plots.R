################ Delta Plots ###################

### Load in Data ########
founderA1_read=read.table("temp/genome/founderA1_10kb_cov20_maf5.txt",header=T)  
founderA2_read=read.table("temp/genome/founderA2_10kb_cov20_maf5.txt",header=T) 
founderB3_read=read.table("temp/genome/founderB3_10kb_cov20_maf5.txt",header=T)  
founderB4_read=read.table("temp/genome/founderB4_10kb_cov20_maf5.txt",header=T) 

### Convert POS to MB for easier visualization with haplotypes 

non_cumulative_offsets <- c(
  230218, 813184, 316620, 1531933, 576874, 270161, 
  1090940, 562643, 439888, 745751, 666816, 1078177, 
  924431, 784333, 1091291
)

# Calculate cumulative offsets
cumulative_offsets <- cumsum(non_cumulative_offsets)

# Add a leading zero for the first chromosome
cumulative_offsets <- c(0, cumulative_offsets)

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
  df$MB_h <- df$pos + offsets[levels(df$chr)][df$chr_numeric]
  
  # Drop the temporary numeric column
  df$chr_numeric <- NULL
  
  return(df)
}

# Apply the function to each dataset with cumulative offsets
founderA1 <- apply_specific_offsets(founderA1_read, named_cumulative_offsets)
founderA2 <- apply_specific_offsets(founderA2_read, named_cumulative_offsets)
founderB3 <- apply_specific_offsets(founderB3_read, named_cumulative_offsets)
founderB4 <- apply_specific_offsets(founderB4_read, named_cumulative_offsets)

################# Untransformed Haplotype Frequencies ########################

#find mean haplotype frequency for young vs old for each position in each FOUNDER. #this finds means across old vs young for each founder, in the code above we found differences for each replicate and then took the mean
founderA1$old.mean=apply(founderA1[4:15],1,mean)
founderA1$young.mean=apply(founderA1[16:27],1,mean)

founderA2$old.mean=apply(founderA2[4:15],1,mean) 
founderA2$young.mean=apply(founderA2[16:27],1,mean)

founderB3$old.mean=apply(founderB3[4:15],1,mean) 
founderB3$young.mean=apply(founderB3[16:27],1,mean) 

founderB4$old.mean=apply(founderB4[4:15],1,mean) 
founderB4$young.mean=apply(founderB4[16:27],1,mean) 

founder_list<-list(founderA1, founderA2, founderB3, founderB4) ### includes all data, but we really only need the means we just generated

results<-list()
for (df in founder_list) {
  sq_diff<- (df$old.mean-df$young.mean)^2
  results[[length(results)+1]]<-sq_diff ##find (h_o,j-h_y,j)^2) from EQ2 for each of four founders
} ### There should be one number for each POS in each haplotype which represents the difference between old and young mean squared 

all_sq_diff<-do.call(cbind, results)

avg_sq_diff<-rowMeans(all_sq_diff) #now we take the sum of all 4 haplotypes at each POS, and divide it my the number of haplotypes (n) to get average
avg_sq_diff_f<-(sqrt(avg_sq_diff)) # now we take the square root (and can multiply by 100 if we want a % value)
MB_h <- founderA1[, "MB_h"] # add the MB scale 
MB_h<-MB_h/1e6 ### convert to correct units

y6 <- cbind(avg_sq_diff_f, MB_h)
y6<-as.data.frame(y6)


tiff(file = "figures/G3_submission/SuppFig8.tiff", height = 4, width = 7,
     units = "in",
     res = 600)

# pdf("figures/Supp_Fig4_Delta.pdf",height=6,width=12)
NT<-plot(y6$MB_h,y6$avg_sq_diff_f,
     type="l",
     ylim=c(0,.1),
     col="black",
     xlab="MB",
     xaxt="n",
     ylab="Delta")

bottom3<-0
top3<-.1

rect(0.23,bottom3,1.04,top3,col="grey80",lty=0)  
rect(1.36,bottom3,2.89,top3,col="grey80",lty=0)  
rect(3.46,bottom3,3.73,top3,col="grey80",lty=0)
rect(4.82,bottom3,5.39,top3,col="grey80",lty=0)  
rect(5.83,bottom3,6.57,top3,col="grey80",lty=0)    
rect(7.24,bottom3,8.32,top3,col="grey80",lty=0)
rect(9.24,bottom3,10.03,top3,col="grey80",lty=0)
rect(11.12,bottom3,12.07,top3,col="grey80",lty=0) 

mtext("C2",line = .5,side=1, at =0.635, cex=1.2)
mtext("C4",line = .5,side=1, at =2.125, cex=1.2)
mtext("C6",line = .5,side=1, at =3.595, cex=1.2)
mtext("C8",line = .5,side=1, at =5.105, cex=1.2)
mtext("C10",line = .5,side=1, at =6.2, cex=1.2)
mtext("C12",line = .5,side=1, at =7.78, cex=1.2)
mtext("C14",line = .5,side=1, at =9.635, cex=1.2)
mtext("C16",line = .5,side=1, at =11.595, cex=1.2)

lines(y6$MB_h,y6$avg_sq_diff_f,col="black")

dev.off()

############## Transformed Haplotype Frequencies ########################################
# 
# #Arcsine sqrt transformation of haplotype frequencies...better approximation of normal distribution
# founderA1_t<-founderA1_read
# founderA2_t<-founderA2_read
# founderB3_t<-founderB3_read
# founderB4_t<-founderB4_read
# 
# # asin sqrt transform all frequencies 
# founderA1_t[, 3:27] <- apply(founderA1_t[, 3:27], 2, function(x) asin(sqrt(x)))
# founderA2_t[, 3:27] <- apply(founderA2_t[, 3:27], 2, function(x) asin(sqrt(x)))
# founderB3_t[, 3:27] <- apply(founderB3_t[, 3:27], 2, function(x) asin(sqrt(x)))
# founderB4_t[, 3:27] <- apply(founderB4_t[, 3:27], 2, function(x) asin(sqrt(x)))
# 
# ### Plot along MB scale #ads MB_h column calculated from above
# founderA1_t <- apply_specific_offsets(founderA1_t, named_cumulative_offsets)
# founderA2_t <- apply_specific_offsets(founderA2_t, named_cumulative_offsets)
# founderB3_t <- apply_specific_offsets(founderB3_t, named_cumulative_offsets)
# founderB4_t <- apply_specific_offsets(founderB4_t, named_cumulative_offsets)
# 
# ####find mean haplotype frequency for young vs old for each position in each FOUNDER. #this finds means across old vs young for each founder, in the code above we found differences for each replicate and then took the mean
# founderA1_t$old.mean=apply(founderA1_t[4:15],1,mean)
# founderA1_t$young.mean=apply(founderA1_t[16:27],1,mean)
# 
# founderA2_t$old.mean=apply(founderA2_t[4:15],1,mean) 
# founderA2_t$young.mean=apply(founderA2_t[16:27],1,mean)
# 
# founderB3_t$old.mean=apply(founderB3_t[4:15],1,mean) 
# founderB3_t$young.mean=apply(founderB3_t[16:27],1,mean) 
# 
# founderB4_t$old.mean=apply(founderB4_t[4:15],1,mean) 
# founderB4_t$young.mean=apply(founderB4_t[16:27],1,mean) 
# 
# founder_list_t<-list(founderA1_t, founderA2_t, founderB3_t, founderB4_t)
# 
# #finds the absolute value of the difference between the old and young transformed mean
# results_t<-list()
# for (df in founder_list_t) {
#   sq_diff<- (df$old.mean-df$young.mean)^2
#   results_t[[length(results_t)+1]]<-sq_diff ##find (h_o,j-h_y,j)^2) from EQ2 for each of four founders
# }
# 
# all_sq_diff_t<-do.call(cbind, results_t)
# avg_sq_diff_t<-rowMeans(all_sq_diff_t) #takes the sum of the four founders and divides by "n" to get average
# delta_t<-(sqrt(avg_sq_diff_t))
# MB_h<- founderA1_t[, "MB_h"]
# MB_h<- MB_h/1e6
# 
# y5 <- cbind(delta_t, MB_h)
# y5<-as.data.frame(y5)
# head(y5)
# 
# pdf("WG_D_transformation.pdf",height=6,width=12)
# TF<-plot(y5$MB_h,y5$delta_t,
#      type="l",
#      ylim=c(0,.2),
#      col="black",
#      xlab="MB",
#      xaxt="n",
#      ylab="diff in hap freq",
#      main="mean diffs hap freqs old-young [1-12]")
# 
# bottom3<-0
# top3<-.2
# 
# rect(0.23,bottom3,1.04,top3,col="grey80",lty=0)  
# rect(1.36,bottom3,2.89,top3,col="grey80",lty=0)  
# rect(3.46,bottom3,3.73,top3,col="grey80",lty=0)
# rect(4.82,bottom3,5.39,top3,col="grey80",lty=0)  
# rect(5.83,bottom3,6.57,top3,col="grey80",lty=0)    
# rect(7.24,bottom3,8.32,top3,col="grey80",lty=0)
# rect(9.24,bottom3,10.03,top3,col="grey80",lty=0)
# rect(11.12,bottom3,12.07,top3,col="grey80",lty=0) 
# 
# mtext("C2",line = .5,side=1, at =0.635, cex=1.2)
# mtext("C4",line = .5,side=1, at =2.125, cex=1.2)
# mtext("C6",line = .5,side=1, at =3.595, cex=1.2)
# mtext("C8",line = .5,side=1, at =5.105, cex=1.2)
# mtext("C10",line = .5,side=1, at =6.2, cex=1.2)
# mtext("C12",line = .5,side=1, at =7.78, cex=1.2)
# mtext("C14",line = .5,side=1, at =9.635, cex=1.2)
# mtext("C16",line = .5,side=1, at =11.595, cex=1.2)
# 
# lines(y5$MB_h,y5$delta_t,col="black")
# 
# dev.off()
# 
# pdf("D_comparison.pdf", width=8, height=8)
# par(mfrow=c(2,1))
# plot(y6$MB_h,y6$avg_sq_diff_f,
#          type="l",
#          ylim=c(0,.1),
#          col="black",
#          xlab="MB",
#          xaxt="n",
#          ylab="Delta",
#          main="mean diffs hap freqs old-young [1-12]")
# 
# bottom3<-0
# top3<-.1
# 
# rect(0.23,bottom3,1.04,top3,col="grey80",lty=0)  
# rect(1.36,bottom3,2.89,top3,col="grey80",lty=0)  
# rect(3.46,bottom3,3.73,top3,col="grey80",lty=0)
# rect(4.82,bottom3,5.39,top3,col="grey80",lty=0)  
# rect(5.83,bottom3,6.57,top3,col="grey80",lty=0)    
# rect(7.24,bottom3,8.32,top3,col="grey80",lty=0)
# rect(9.24,bottom3,10.03,top3,col="grey80",lty=0)
# rect(11.12,bottom3,12.07,top3,col="grey80",lty=0) 
# 
# mtext("C2",line = .5,side=1, at =0.635, cex=1.2)
# mtext("C4",line = .5,side=1, at =2.125, cex=1.2)
# mtext("C6",line = .5,side=1, at =3.595, cex=1.2)
# mtext("C8",line = .5,side=1, at =5.105, cex=1.2)
# mtext("C10",line = .5,side=1, at =6.2, cex=1.2)
# mtext("C12",line = .5,side=1, at =7.78, cex=1.2)
# mtext("C14",line = .5,side=1, at =9.635, cex=1.2)
# mtext("C16",line = .5,side=1, at =11.595, cex=1.2)
# 
# lines(y6$MB_h,y6$avg_sq_diff_f,col="black")
# 
# plot(y5$MB_h,y5$delta_t,
#          type="l",
#          ylim=c(0,.16),
#          col="black",
#          xlab="MB",
#          xaxt="n",
#          ylab="Transformed Delta",
#          main="mean diffs hap freqs old-young [1-12]")
# 
# bottom3<-0
# top3<-.16
# 
# rect(0.23,bottom3,1.04,top3,col="grey80",lty=0)  
# rect(1.36,bottom3,2.89,top3,col="grey80",lty=0)  
# rect(3.46,bottom3,3.73,top3,col="grey80",lty=0)
# rect(4.82,bottom3,5.39,top3,col="grey80",lty=0)  
# rect(5.83,bottom3,6.57,top3,col="grey80",lty=0)    
# rect(7.24,bottom3,8.32,top3,col="grey80",lty=0)
# rect(9.24,bottom3,10.03,top3,col="grey80",lty=0)
# rect(11.12,bottom3,12.07,top3,col="grey80",lty=0) 
# 
# mtext("C2",line = .5,side=1, at =0.635, cex=1.2)
# mtext("C4",line = .5,side=1, at =2.125, cex=1.2)
# mtext("C6",line = .5,side=1, at =3.595, cex=1.2)
# mtext("C8",line = .5,side=1, at =5.105, cex=1.2)
# mtext("C10",line = .5,side=1, at =6.2, cex=1.2)
# mtext("C12",line = .5,side=1, at =7.78, cex=1.2)
# mtext("C14",line = .5,side=1, at =9.635, cex=1.2)
# mtext("C16",line = .5,side=1, at =11.595, cex=1.2)
# 
# lines(y5$MB_h,y5$delta_t,col="black")
# 
# dev.off()
# 
# ################### Look at haplotypes too ##########################
# 
# for (i in 1:12) {
#   founderA1[[paste0("diff", i)]] <- founderA1[[paste0("hap.old", sprintf("%02d", i))]] - founderA1[[paste0("hap.young", sprintf("%02d", i))]]
# }
# 
# for (i in 1:12) {
#   founderA2[[paste0("diff", i)]] <- founderA2[[paste0("hap.old", sprintf("%02d", i))]] - founderA2[[paste0("hap.young", sprintf("%02d", i))]]
# }
# 
# for (i in 1:12) {
#   founderB3[[paste0("diff", i)]] <- founderB3[[paste0("hap.old", sprintf("%02d", i))]] - founderB3[[paste0("hap.young", sprintf("%02d", i))]]
# }
# 
# for (i in 1:12) {
#   founderB4[[paste0("diff", i)]] <- founderB4[[paste0("hap.old", sprintf("%02d", i))]] - founderB4[[paste0("hap.young", sprintf("%02d", i))]]
# }
# 
# View(founderA1[,c(31:42)])
# ## Find the mean differences across all pairs from the differences calculated above, and the variance
# founderA1$meandiff=apply(founderA1[31:42],1,mean) # this is a first-pass simple effort to quantify consistent differences in haplotypes across replicates, per fou31:42
# founderA2$meandiff=apply(founderA2[31:42],1,mean) # this is a first-pass simple effort to quantify consistent differences in haplotypes across replicates, per fou31:42
# founderB3$meandiff=apply(founderB3[31:42],1,mean) # this is a first-pass simple effort to quantify consistent differences in haplotypes across replicates, per fou31:42
# founderB4$meandiff=apply(founderB4[31:42],1,mean) # this is a first-pass simple effort to quantify consistent differences in haplotypes across replicates, per founder
# 
# top6<-0.1
# bottom<-0
#   
# View(founderA1)
# plot(founderA1$MB_h,founderA1$meandiff,
#      type="l",
#      xlab= "Genomic Position",
#      ylim=c(bottom,top6),
#      col="darkgreen",
#      xaxt= "n",
#      yaxt= "n",
#      cex.lab= 2,
#      ylab=bquote("Diff in Hap Freq"))
# 
# mtext("B", side = 3, adj = -0.1, line = 1.5, cex = 2.2, font=2)
# 
# y_ticks <- seq(0, top6, length.out = 5)  # Modify as needed
# y_labels <- round(y_ticks, 2)  # Adjust the number of decimals for your labels
# 
# 
# rect(0.23,bottom,1.04,top6,col="grey80",lty=0)  
# rect(1.36,bottom,2.89,top6,col="grey80",lty=0)  
# rect(3.46,bottom,3.73,top6,col="grey80",lty=0)
# rect(4.82,bottom,5.39,top6,col="grey80",lty=0)  
# rect(5.83,bottom,6.57,top6,col="grey80",lty=0)    
# rect(7.24,bottom,8.32,top6,col="grey80",lty=0)
# rect(9.24,bottom,10.03,top6,col="grey80",lty=0)
# rect(11.12,bottom,12.07,top6,col="grey80",lty=0) 
# 
# lines(founderA1$MB,founderA1$meandiff,col="darkgreen")
# lines(founderA2$MB,founderA2$meandiff,col="red")
# lines(founderB3$MB,founderB3$meandiff,col="goldenrod")
# lines(founderB4$MB,founderB4$meandiff,col="blue") # these colors have been used before to refer to these specific founders
# 
# #add chromosome labels at midpts
# mtext("C2",line = .5,side=1, at =0.635, cex=1.3)
# mtext("C4",line = .5,side=1, at =2.125, cex=1.3)
# mtext("C6",line = .5,side=1, at =3.595, cex=1.3)
# mtext("C8",line = .5,side=1, at =5.105, cex=1.3)
# mtext("C10",line = .5,side=1, at =6.2, cex=1.3)
# mtext("C12",line = .5,side=1, at =7.78, cex=1.3)
# mtext("C14",line = .5,side=1, at =9.635, cex=1.3)
# mtext("C16",line = .5,side=1, at =11.595, cex=1.3)
# 
# legend("topright", 
#        legend = c("DBVPG6765", "DBPVG6044","YPS128", "Y12"),    # Names of the data
#        col = c("darkgreen", "red", "goldenrod", "blue"),          # Corresponding colors
#        lty = c(1, 1, 1, 1),                 # Corresponding point shapes
#        cex = 1.1, 
#        lwd= 2,
#        bg="white", 
#        ncol=2)    
# 
# axis(2, at = y_ticks, labels = y_labels, cex.axis = 1.8)
# 
# 
# #############
# 
# 
# pdf("D_comparison_haps.pdf", width=8, height=12)
# par(mfrow=c(3,1))
# plot(y6$MB_h,y6$avg_sq_diff_f,
#      type="l",
#      ylim=c(0,.1),
#      col="black",
#      xlab="MB",
#      xaxt="n",
#      ylab="Delta",
#      main="mean diffs hap freqs old-young [1-12]")
# 
# bottom3<-0
# top3<-.1
# 
# rect(0.23,bottom3,1.04,top3,col="grey80",lty=0)  
# rect(1.36,bottom3,2.89,top3,col="grey80",lty=0)  
# rect(3.46,bottom3,3.73,top3,col="grey80",lty=0)
# rect(4.82,bottom3,5.39,top3,col="grey80",lty=0)  
# rect(5.83,bottom3,6.57,top3,col="grey80",lty=0)    
# rect(7.24,bottom3,8.32,top3,col="grey80",lty=0)
# rect(9.24,bottom3,10.03,top3,col="grey80",lty=0)
# rect(11.12,bottom3,12.07,top3,col="grey80",lty=0) 
# 
# mtext("C2",line = .5,side=1, at =0.635, cex=1.2)
# mtext("C4",line = .5,side=1, at =2.125, cex=1.2)
# mtext("C6",line = .5,side=1, at =3.595, cex=1.2)
# mtext("C8",line = .5,side=1, at =5.105, cex=1.2)
# mtext("C10",line = .5,side=1, at =6.2, cex=1.2)
# mtext("C12",line = .5,side=1, at =7.78, cex=1.2)
# mtext("C14",line = .5,side=1, at =9.635, cex=1.2)
# mtext("C16",line = .5,side=1, at =11.595, cex=1.2)
# 
# lines(y6$MB_h,y6$avg_sq_diff_f,col="black")
# 
# plot(founderA1$MB_h,founderA1$meandiff,
#      type="l",
#      xlab= "Genomic Position",
#      ylim=c(bottom,top6),
#      col="darkgreen",
#      xaxt= "n",
#      yaxt= "n",
#      cex.lab= 2,
#      ylab=bquote("Diff in Hap Freq"))
# 
# mtext("B", side = 3, adj = -0.1, line = 1.5, cex = 2.2, font=2)
# 
# y_ticks <- seq(0, top6, length.out = 5)  # Modify as needed
# y_labels <- round(y_ticks, 2)  # Adjust the number of decimals for your labels
# 
# 
# rect(0.23,bottom,1.04,top6,col="grey80",lty=0)  
# rect(1.36,bottom,2.89,top6,col="grey80",lty=0)  
# rect(3.46,bottom,3.73,top6,col="grey80",lty=0)
# rect(4.82,bottom,5.39,top6,col="grey80",lty=0)  
# rect(5.83,bottom,6.57,top6,col="grey80",lty=0)    
# rect(7.24,bottom,8.32,top6,col="grey80",lty=0)
# rect(9.24,bottom,10.03,top6,col="grey80",lty=0)
# rect(11.12,bottom,12.07,top6,col="grey80",lty=0) 
# 
# lines(founderA1$MB,founderA1$meandiff,col="darkgreen")
# lines(founderA2$MB,founderA2$meandiff,col="red")
# lines(founderB3$MB,founderB3$meandiff,col="goldenrod")
# lines(founderB4$MB,founderB4$meandiff,col="blue") # these colors have been used before to refer to these specific founders
# 
# #add chromosome labels at midpts
# mtext("C2",line = .5,side=1, at =0.635, cex=1.3)
# mtext("C4",line = .5,side=1, at =2.125, cex=1.3)
# mtext("C6",line = .5,side=1, at =3.595, cex=1.3)
# mtext("C8",line = .5,side=1, at =5.105, cex=1.3)
# mtext("C10",line = .5,side=1, at =6.2, cex=1.3)
# mtext("C12",line = .5,side=1, at =7.78, cex=1.3)
# mtext("C14",line = .5,side=1, at =9.635, cex=1.3)
# mtext("C16",line = .5,side=1, at =11.595, cex=1.3)
# 
# legend("topright", 
#        legend = c("DBVPG6765", "DBPVG6044","YPS128", "Y12"),    # Names of the data
#        col = c("darkgreen", "red", "goldenrod", "blue"),          # Corresponding colors
#        lty = c(1, 1, 1, 1),                 # Corresponding point shapes
#        cex = 1.1, 
#        lwd= 2,
#        bg="white", 
#        ncol=2)    
# 
# axis(2, at = y_ticks, labels = y_labels, cex.axis = 1.8)
# 
# plot(y5$MB_h,y5$delta_t,
#      type="l",
#      ylim=c(0,.16),
#      col="black",
#      xlab="MB",
#      xaxt="n",
#      ylab="Transformed Delta",
#      main="mean diffs hap freqs old-young [1-12]")
# 
# bottom3<-0
# top3<-.16
# 
# rect(0.23,bottom3,1.04,top3,col="grey80",lty=0)  
# rect(1.36,bottom3,2.89,top3,col="grey80",lty=0)  
# rect(3.46,bottom3,3.73,top3,col="grey80",lty=0)
# rect(4.82,bottom3,5.39,top3,col="grey80",lty=0)  
# rect(5.83,bottom3,6.57,top3,col="grey80",lty=0)    
# rect(7.24,bottom3,8.32,top3,col="grey80",lty=0)
# rect(9.24,bottom3,10.03,top3,col="grey80",lty=0)
# rect(11.12,bottom3,12.07,top3,col="grey80",lty=0) 
# 
# mtext("C2",line = .5,side=1, at =0.635, cex=1.2)
# mtext("C4",line = .5,side=1, at =2.125, cex=1.2)
# mtext("C6",line = .5,side=1, at =3.595, cex=1.2)
# mtext("C8",line = .5,side=1, at =5.105, cex=1.2)
# mtext("C10",line = .5,side=1, at =6.2, cex=1.2)
# mtext("C12",line = .5,side=1, at =7.78, cex=1.2)
# mtext("C14",line = .5,side=1, at =9.635, cex=1.2)
# mtext("C16",line = .5,side=1, at =11.595, cex=1.2)
# 
# lines(y5$MB_h,y5$delta_t,col="black")
# 
# dev.off()