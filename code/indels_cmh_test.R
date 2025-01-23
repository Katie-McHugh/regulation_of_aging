### Indels CMH test

indels3<-read.csv("temp/indels_cov20_maf05.csv")

#create a bunch of vectors
old_alt_01<-indels3$alt_old_01
old_alt_02<-indels3$alt_old_02
old_alt_03<-indels3$alt_old_03
old_alt_04<-indels3$alt_old_04
old_alt_05<-indels3$alt_old_05
old_alt_06<-indels3$alt_old_06
old_alt_07<-indels3$alt_old_07
old_alt_08<-indels3$alt_old_08
old_alt_09<-indels3$alt_old_09
old_alt_10<-indels3$alt_old_10
old_alt_11<-indels3$alt_old_11
old_alt_12<-indels3$alt_old_12
old_ref_01<-indels3$N_old_01-indels3$alt_old_01
old_ref_02<-indels3$N_old_02-indels3$alt_old_02
old_ref_03<-indels3$N_old_03-indels3$alt_old_03
old_ref_04<-indels3$N_old_04-indels3$alt_old_04
old_ref_05<-indels3$N_old_05-indels3$alt_old_05
old_ref_06<-indels3$N_old_06-indels3$alt_old_06
old_ref_07<-indels3$N_old_07-indels3$alt_old_07
old_ref_08<-indels3$N_old_08-indels3$alt_old_08
old_ref_09<-indels3$N_old_09-indels3$alt_old_09
old_ref_10<-indels3$N_old_10-indels3$alt_old_10
old_ref_11<-indels3$N_old_11-indels3$alt_old_11
old_ref_12<-indels3$N_old_12-indels3$alt_old_12

young_alt_01<-indels3$alt_young_01
young_alt_02<-indels3$alt_young_02
young_alt_03<-indels3$alt_young_03
young_alt_04<-indels3$alt_young_04
young_alt_05<-indels3$alt_young_05
young_alt_06<-indels3$alt_young_06
young_alt_07<-indels3$alt_young_07
young_alt_08<-indels3$alt_young_08
young_alt_09<-indels3$alt_young_09
young_alt_10<-indels3$alt_young_10
young_alt_11<-indels3$alt_young_11
young_alt_12<-indels3$alt_young_12
young_ref_01<-indels3$N_young_01-indels3$alt_young_01
young_ref_02<-indels3$N_young_02-indels3$alt_young_02
young_ref_03<-indels3$N_young_03-indels3$alt_young_03
young_ref_04<-indels3$N_young_04-indels3$alt_young_04
young_ref_05<-indels3$N_young_05-indels3$alt_young_05
young_ref_06<-indels3$N_young_06-indels3$alt_young_06
young_ref_07<-indels3$N_young_07-indels3$alt_young_07
young_ref_08<-indels3$N_young_08-indels3$alt_young_08
young_ref_09<-indels3$N_young_09-indels3$alt_young_09
young_ref_10<-indels3$N_young_10-indels3$alt_young_10
young_ref_11<-indels3$N_young_11-indels3$alt_young_11
young_ref_12<-indels3$N_young_12-indels3$alt_young_12

### use all 12
#create my matrices
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
indels3$logp<-logp
indels3$pval<-pval

# Save the file to the new directory
write.csv(indels3, file = "temp/indels_CMHtest_results.csv", row.names = FALSE)
