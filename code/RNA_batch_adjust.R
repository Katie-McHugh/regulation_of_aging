### Correct for batch effect from 2 different sorting groups

#Read in Gene Count Matrix, subset
counts <- read.table("data/gene_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
CountData<-counts[,-(1:5)]

## read in design file and designate as factors/matrix
colData <- read.table("data/design.txt", header = TRUE, sep = "\t", row.names = 1)
colData<-as.matrix(colData)
View(colData)

# batch information
batch <- c(rep(1, 6), rep(2, 6), rep(1, 6), rep(2, 6)) #Organize data into batch 1 (pairs 1-6, sorted on Day 1) and batch 2 (pairs 7-12, sorted on Day 2)
age<-c(rep(1, 12), rep(2, 12)) # indicate whether the replicate is "old" (1) or "young" (2)
#pair<-rep(1:12, 2) #indicate the paired nature of the old and young replicates # don't include this in CombatSeq

# Run ComBat_seq
adj_counts <- ComBat_seq(CountData, batch=batch, group=age) #ignore the covariates...package is not super clear on what they are used for, and we don't want to include pair in both ComBat-seq and the DEseq model

write.csv(adj_counts, file="temp/gcm_combatseq.csv")
