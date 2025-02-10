### Run a CMH test on SNP data

snps3<-read.csv("temp/genome/GWAS_SNPS_cov20_maf5.csv")

old_alt_01<-snps3$alt_old_01
old_alt_02<-snps3$alt_old_02
old_alt_03<-snps3$alt_old_03
old_alt_04<-snps3$alt_old_04
old_alt_05<-snps3$alt_old_05
old_alt_06<-snps3$alt_old_06
old_alt_07<-snps3$alt_old_07
old_alt_08<-snps3$alt_old_08
old_alt_09<-snps3$alt_old_09
old_alt_10<-snps3$alt_old_10
old_alt_11<-snps3$alt_old_11
old_alt_12<-snps3$alt_old_12
old_ref_01<-snps3$N_old_01-snps3$alt_old_01
old_ref_02<-snps3$N_old_02-snps3$alt_old_02
old_ref_03<-snps3$N_old_03-snps3$alt_old_03
old_ref_04<-snps3$N_old_04-snps3$alt_old_04
old_ref_05<-snps3$N_old_05-snps3$alt_old_05
old_ref_06<-snps3$N_old_06-snps3$alt_old_06
old_ref_07<-snps3$N_old_07-snps3$alt_old_07
old_ref_08<-snps3$N_old_08-snps3$alt_old_08
old_ref_09<-snps3$N_old_09-snps3$alt_old_09
old_ref_10<-snps3$N_old_10-snps3$alt_old_10
old_ref_11<-snps3$N_old_11-snps3$alt_old_11
old_ref_12<-snps3$N_old_12-snps3$alt_old_12
young_alt_01<-snps3$alt_young_01
young_alt_02<-snps3$alt_young_02
young_alt_03<-snps3$alt_young_03
young_alt_04<-snps3$alt_young_04
young_alt_05<-snps3$alt_young_05
young_alt_06<-snps3$alt_young_06
young_alt_07<-snps3$alt_young_07
young_alt_08<-snps3$alt_young_08
young_alt_09<-snps3$alt_young_09
young_alt_10<-snps3$alt_young_10
young_alt_11<-snps3$alt_young_11
young_alt_12<-snps3$alt_young_12
young_ref_01<-snps3$N_young_01-snps3$alt_young_01
young_ref_02<-snps3$N_young_02-snps3$alt_young_02
young_ref_03<-snps3$N_young_03-snps3$alt_young_03
young_ref_04<-snps3$N_young_04-snps3$alt_young_04
young_ref_05<-snps3$N_young_05-snps3$alt_young_05
young_ref_06<-snps3$N_young_06-snps3$alt_young_06
young_ref_07<-snps3$N_young_07-snps3$alt_young_07
young_ref_08<-snps3$N_young_08-snps3$alt_young_08
young_ref_09<-snps3$N_young_09-snps3$alt_young_09
young_ref_10<-snps3$N_young_10-snps3$alt_young_10
young_ref_11<-snps3$N_young_11-snps3$alt_young_11
young_ref_12<-snps3$N_young_12-snps3$alt_young_12

### use all 12
A0_12<-rbind(old_alt_01, old_alt_02, old_alt_03, old_alt_04, old_alt_05, old_alt_06, old_alt_07, old_alt_08, old_alt_09, old_alt_10, old_alt_11, old_alt_12 )

a0_12<-rbind(old_ref_01, old_ref_02, old_ref_03, old_ref_04, old_ref_05, old_ref_06, old_ref_07, old_ref_08, old_ref_09, old_ref_10, old_ref_11, old_ref_12)

At_12<-rbind(young_alt_01, young_alt_02, young_alt_03, young_alt_04, young_alt_05, young_alt_06, young_alt_07, young_alt_08, young_alt_09, young_alt_10, young_alt_11, young_alt_12)

at_12<-rbind(young_ref_01, young_ref_02, young_ref_03, young_ref_04, young_ref_05, young_ref_06, young_ref_07, young_ref_08, young_ref_09, young_ref_10, young_ref_11, young_ref_12)


deal_with_nas<-0.01
A3<-A0_12+deal_with_nas
a3<-a0_12+deal_with_nas
A4<-At_12+deal_with_nas
a4<-at_12+deal_with_nas

logp<-cmh.test(A3, a3, A4, a4, log="TRUE") 
pval<-cmh.test(A3, a3, A4, a4, log="FALSE")
snps3$logp<-logp
snps3$pval<-pval


#### Prepare Genome Scale
snps4<-subset(snps3, CHROM!="chrmito")
## GW plotting code (old school plot, can be converted to ggplot.  This corresponds to the nuclear genome scale for the R64 genome release)
Gaxis <- numeric(length=0)
chrs<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16')
for (k in 1:length(chrs)){
  data.samp<-subset(snps4,CHROM==chrs[k])
  if (chrs[k]=='chr1'){Gaxis.samp<-data.samp$POS}
  if (chrs[k]=='chr2'){Gaxis.samp<-data.samp$POS+230218}
  if (chrs[k]=='chr3'){Gaxis.samp<-data.samp$POS+230218+813184}
  if (chrs[k]=='chr4'){Gaxis.samp<-data.samp$POS+230218+813184+316620}
  if (chrs[k]=='chr5'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933}
  if (chrs[k]=='chr6'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874}
  if (chrs[k]=='chr7'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='chr8'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='chr9'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='chr10'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='chr11'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='chr12'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='chr13'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='chr14'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='chr15'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='chr16'){Gaxis.samp<-data.samp$POS+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}

MB=Gaxis/1e6  ## MB stands for megabases, we are dividing by 1 million to put things on a megabase scale
snps4$MB=MB ## now we have a vector in our SNP table that we can use as a x-axis variable for plotting

# Save the file to the new directory
write.csv(snps3, file = "temp/genome/WG_CMHtest_results.csv", row.names = FALSE)
write.csv(snps4, file = "temp/genome/WG_CMHtest_results_nuclear.csv", row.names = FALSE)
