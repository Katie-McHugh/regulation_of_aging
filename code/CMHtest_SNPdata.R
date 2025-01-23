### Run a CMH test on SNP data

snps3<-read.csv("temp/GWAS_SNPS_cov20_maf5.csv")

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

# Specify the file path in the new directory
file_path2 <- file.path("temp/WG_CMHtest_results.csv")

# Save the file to the new directory
write.csv(snps3, file = file_path2, row.names = FALSE)