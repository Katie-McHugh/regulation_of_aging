### Haplotype estimation full script. ## this was all done before moving the script 
# to the respository. It was run and saved, this was the code used.

# load filtered data from Analysis_eNotebook.Rmd file
snps3=read.csv("temp/genome/GWAS_SNPS_cov20_maf5.csv", header=TRUE)

## haplotype estimator
founders=read.table("data/founder_states.txt",header=T)
merge=merge(founders,snps3) #will use current filtering parameters

freq <- apply(merge[,10:62],1,function(x) (x[seq(1,48,2)]/(x[seq(2,49,2)]))) 
freq2 <- t(freq)  # this transposes the new matrix
index=c(1,2,5:9)
freq3 <- cbind(merge[index],freq2)  

## this is the table we will use to run the haplotype estimator.  You could save this for the future so that you don't have to read it back in - but I wanted to make clear notes about how I generated this table (as we may want to tweak it in the future)
write.table(freq3,file="temp/genome/GWAS_snps_with_founders.txt",quote=FALSE)


#######     HAPLOTYPE ESTIMATOR   #####################################################################################################
#####################################################################################################################

maxfun <- function(p,X,Y,N){
  SS1 <- sum((Y - X %*% p)^2)
  w <- 100
  SS2 <- w * N * ((sum(p) - 1)^2)
  SS1 + SS2
}

gethaps <- function(x,mafs,pos,fstates){
  half.window.size <<- 5000
  window <- (pos > (x - half.window.size)) & (pos < (x + half.window.size))
  N <<- sum(window)
  Y <<- as.matrix(mafs[window],ncol=1)
  X <<- as.matrix(fstates[window,],ncol=4)
  p <<- as.matrix(rep(1/4,4),ncol=1)
  out <- optim(p, maxfun, X=X,Y=Y,N=N, method = "L-BFGS-B", lower = rep(0.001,4))
  # constraining with an upper bound causes problems, excess of p=0 values
  # 4 haplotype frequencies
  out$par
}

##########################################################################

haps=freq3

## you have to do the haplotype estimation chromosome by chromosome
## we also have to choose a window size.  starting with 10kb (half window = 5000 above) but this can be messed with


#### CHROMOSOME 1

xx=subset(haps,CHROM=="chr1")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr1",length(testpositions))

pos=testpositions

chr1_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr1_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr1_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr1_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr1_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr1_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr1_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr1_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

############### CHROMOSOME 2 ########################

xx=subset(haps,CHROM=="chr2")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr2",length(testpositions))

pos=testpositions

chr2_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr2_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr2_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr2_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr2_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr2_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr2_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr2_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

##################### CHROMOSOME 3 ##################

xx=subset(haps,CHROM=="chr3")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr3",length(testpositions))

pos=testpositions

chr3_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr3_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr3_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr3_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr3_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr3_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr3_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr3_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

################### CHROMOSOME 4 #########################

xx=subset(haps,CHROM=="chr4")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr4",length(testpositions))

pos=testpositions

chr4_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr4_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr4_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr4_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr4_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr4_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr4_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr4_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

########################### CHROMOSOME 5 #######################################

xx=subset(haps,CHROM=="chr5")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr5",length(testpositions))

pos=testpositions

chr5_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr5_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr5_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr5_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr5_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr5_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr5_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr5_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

######################### CHROMOSOME 6 ###################

xx=subset(haps,CHROM=="chr6")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr6",length(testpositions))

pos=testpositions

chr6_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr6_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr6_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr6_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr6_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr6_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr6_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr6_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

########## CHROMOSOME 7 ####################

xx=subset(haps,CHROM=="chr7")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr7",length(testpositions))

pos=testpositions

chr7_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr7_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr7_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr7_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr7_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr7_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr7_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr7_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

############################# CHROMOSOME 8 ############################

xx=subset(haps,CHROM=="chr8")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr8",length(testpositions))

pos=testpositions

chr8_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr8_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr8_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr8_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr8_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr8_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr8_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr8_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

####################### CHROMOSOME 9 #######################

xx=subset(haps,CHROM=="chr9")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr9",length(testpositions))

pos=testpositions

chr9_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr9_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr9_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr9_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr9_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr9_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr9_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr9_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

########### CHROMOSOME 10 ###################


xx=subset(haps,CHROM=="chr10")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr10",length(testpositions))

pos=testpositions

chr10_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr10_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr10_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr10_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr10_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr10_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr10_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr10_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

############## CHROMOSOME 11 ########################

xx=subset(haps,CHROM=="chr11")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr11",length(testpositions))

pos=testpositions

chr11_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr11_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr11_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr11_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr11_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr11_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr11_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr11_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

################# CHROMOSOME 12 ######################

xx=subset(haps,CHROM=="chr12")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr12",length(testpositions))

pos=testpositions

chr12_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr12_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr12_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr12_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr12_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr12_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr12_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr12_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

####################### CHROMOSOME 13 ############################

xx=subset(haps,CHROM=="chr13")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr13",length(testpositions))

pos=testpositions

chr13_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr13_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr13_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr13_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr13_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr13_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr13_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr13_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

############ CHROMOSOME 14 ##########################


xx=subset(haps,CHROM=="chr14")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr14",length(testpositions))

pos=testpositions

chr14_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr14_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr14_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr14_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr14_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr14_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr14_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr14_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

################### CHROMOSOME 15 ########################


xx=subset(haps,CHROM=="chr15")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr15",length(testpositions))

pos=testpositions

chr15_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr15_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr15_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr15_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr15_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr15_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr15_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr15_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

################# CHROMOSOME 16 #########################

xx=subset(haps,CHROM=="chr16")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chr16",length(testpositions))

pos=testpositions

chr16_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chr16_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr16_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chr16_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr16_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chr16_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chr16_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chr16_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

################### MITOCHONDRIA ###################


xx=subset(haps,CHROM=="chrmito")
stepsize=2000 # this can also be messed with
min=min(xx$POS) + 2500 # this is also arbitrary, should be chosen with respect to the step size, just prevents unequal averaging at the chromosome ends
max=max(xx$POS) - 2500

testpositions <- seq(min,max,stepsize)

haps.anc <- t(sapply(testpositions,function(x) gethaps(x,xx[,7],xx$POS,xx[,3:6])))
haps.old.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,8],xx$POS,xx[,3:6])))
haps.old.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,9],xx$POS,xx[,3:6])))
haps.old.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,10],xx$POS,xx[,3:6])))
haps.old.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,11],xx$POS,xx[,3:6])))
haps.old.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,12],xx$POS,xx[,3:6])))
haps.old.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,13],xx$POS,xx[,3:6])))
haps.old.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,14],xx$POS,xx[,3:6])))
haps.old.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,15],xx$POS,xx[,3:6])))
haps.old.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,16],xx$POS,xx[,3:6])))
haps.old.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,17],xx$POS,xx[,3:6])))
haps.old.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,18],xx$POS,xx[,3:6])))
haps.old.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,19],xx$POS,xx[,3:6])))
haps.young.01 <- t(sapply(testpositions,function(x) gethaps(x,xx[,20],xx$POS,xx[,3:6])))
haps.young.02 <- t(sapply(testpositions,function(x) gethaps(x,xx[,21],xx$POS,xx[,3:6])))
haps.young.03 <- t(sapply(testpositions,function(x) gethaps(x,xx[,22],xx$POS,xx[,3:6])))
haps.young.04 <- t(sapply(testpositions,function(x) gethaps(x,xx[,23],xx$POS,xx[,3:6])))
haps.young.05 <- t(sapply(testpositions,function(x) gethaps(x,xx[,24],xx$POS,xx[,3:6])))
haps.young.06 <- t(sapply(testpositions,function(x) gethaps(x,xx[,25],xx$POS,xx[,3:6])))
haps.young.07 <- t(sapply(testpositions,function(x) gethaps(x,xx[,26],xx$POS,xx[,3:6])))
haps.young.08 <- t(sapply(testpositions,function(x) gethaps(x,xx[,27],xx$POS,xx[,3:6])))
haps.young.09 <- t(sapply(testpositions,function(x) gethaps(x,xx[,28],xx$POS,xx[,3:6])))
haps.young.10 <- t(sapply(testpositions,function(x) gethaps(x,xx[,29],xx$POS,xx[,3:6])))
haps.young.11 <- t(sapply(testpositions,function(x) gethaps(x,xx[,30],xx$POS,xx[,3:6])))
haps.young.12 <- t(sapply(testpositions,function(x) gethaps(x,xx[,31],xx$POS,xx[,3:6])))

chr=rep("chrmito",length(testpositions))

pos=testpositions

chrmito_A1=data.frame(chr,pos,haps.anc[,1],haps.old.01[,1],haps.old.02[,1],haps.old.03[,1],haps.old.04[,1],haps.old.05[,1],haps.old.06[,1],haps.old.07[,1],haps.old.08[,1],haps.old.09[,1],haps.old.10[,1],haps.old.11[,1],haps.old.12[,1],haps.young.01[,1],haps.young.02[,1],haps.young.03[,1],haps.young.04[,1],haps.young.05[,1],haps.young.06[,1],haps.young.07[,1],haps.young.08[,1],haps.young.09[,1],haps.young.10[,1],haps.young.11[,1],haps.young.12[,1])
colnames(chrmito_A1)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chrmito_A2=data.frame(chr,pos,haps.anc[,2],haps.old.01[,2],haps.old.02[,2],haps.old.03[,2],haps.old.04[,2],haps.old.05[,2],haps.old.06[,2],haps.old.07[,2],haps.old.08[,2],haps.old.09[,2],haps.old.10[,2],haps.old.11[,2],haps.old.12[,2],haps.young.01[,2],haps.young.02[,2],haps.young.03[,2],haps.young.04[,2],haps.young.05[,2],haps.young.06[,2],haps.young.07[,2],haps.young.08[,2],haps.young.09[,2],haps.young.10[,2],haps.young.11[,2],haps.young.12[,2])
colnames(chrmito_A2)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chrmito_B3=data.frame(chr,pos,haps.anc[,3],haps.old.01[,3],haps.old.02[,3],haps.old.03[,3],haps.old.04[,3],haps.old.05[,3],haps.old.06[,3],haps.old.07[,3],haps.old.08[,3],haps.old.09[,3],haps.old.10[,3],haps.old.11[,3],haps.old.12[,3],haps.young.01[,3],haps.young.02[,3],haps.young.03[,3],haps.young.04[,3],haps.young.05[,3],haps.young.06[,3],haps.young.07[,3],haps.young.08[,3],haps.young.09[,3],haps.young.10[,3],haps.young.11[,3],haps.young.12[,3])
colnames(chrmito_B3)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

chrmito_B4=data.frame(chr,pos,haps.anc[,4],haps.old.01[,4],haps.old.02[,4],haps.old.03[,4],haps.old.04[,4],haps.old.05[,4],haps.old.06[,4],haps.old.07[,4],haps.old.08[,4],haps.old.09[,4],haps.old.10[,4],haps.old.11[,4],haps.old.12[,4],haps.young.01[,4],haps.young.02[,4],haps.young.03[,4],haps.young.04[,4],haps.young.05[,4],haps.young.06[,4],haps.young.07[,4],haps.young.08[,4],haps.young.09[,4],haps.young.10[,4],haps.young.11[,4],haps.young.12[,4])
colnames(chrmito_B4)= c("chr","pos","hap.anc","hap.old01","hap.old02","hap.old03","hap.old04","hap.old05","hap.old06","hap.old07","hap.old08","hap.old09","hap.old10","hap.old11","hap.old12","hap.young01","hap.young02","hap.young03","hap.young04","hap.young05","hap.young06","hap.young07","hap.young08","hap.young09","hap.young10","hap.young11","hap.young12")

############################################################################################################################################################
#####################################################################################################################

## above is how you generate the haplotypes for a single chromosome.  What I did next was to repeat in a loop for all chromosomes, though I did not save that code.  Once I had all chromomosomes saved, I bound them together like this:
founderA1=rbind(chr1_A1,chr2_A1,chr3_A1,chr4_A1,chr5_A1,chr6_A1,chr7_A1,chr8_A1,chr9_A1,chr10_A1,chr11_A1,chr12_A1,chr13_A1,chr14_A1,chr15_A1,chr16_A1)
founderA2=rbind(chr1_A2,chr2_A2,chr3_A2,chr4_A2,chr5_A2,chr6_A2,chr7_A2,chr8_A2,chr9_A2,chr10_A2,chr11_A2,chr12_A2,chr13_A2,chr14_A2,chr15_A2,chr16_A2)
founderB3=rbind(chr1_B3,chr2_B3,chr3_B3,chr4_B3,chr5_B3,chr6_B3,chr7_B3,chr8_B3,chr9_B3,chr10_B3,chr11_B3,chr12_B3,chr13_B3,chr14_B3,chr15_B3,chr16_B3)
founderB4=rbind(chr1_B4,chr2_B4,chr3_B4,chr4_B4,chr5_B4,chr6_B4,chr7_B4,chr8_B4,chr9_B4,chr10_B4,chr11_B4,chr12_B4,chr13_B4,chr14_B4,chr15_B4,chr16_B4)

### Write tables # Would need to reset these
write.table (founderA1,file="founderA1_10kb_cov20_maf5.txt",quote=FALSE)
write.table (founderA2,file="founderA2_10kb_cov20_maf5.txt",quote=FALSE)
write.table (founderB3,file="founderB3_10kb_cov20_maf5.txt",quote=FALSE)
write.table (founderB4,file="founderB4_10kb_cov20_maf5.txt",quote=FALSE)
