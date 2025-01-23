founders_4S=read.table("all_founders.txt",header=T)

index=c(2,3,4,5,10,11,12,13,30,31,32,33,76,77)
test=founders_4S[index]

test$A1=test$alt_YEE_hap_A1_00/test$N_YEE_hap_A1_00
test$A2=test$alt_YEE_hap_A2_00/test$N_YEE_hap_A2_00
test$B3=test$alt_YEE_hap_B3_00/test$N_YEE_hap_B3_00
test$B4=test$alt_YEE_hap_B4_00/test$N_YEE_hap_B4_00
test$anc=test$alt_YEE_rec_4SH_12/test$N_YEE_rec_4SH_12

test2=subset(test,A1==0 | A1==1)
test3=subset(test2,A2==0 | A2==1)
test4=subset(test3,B3==0 | B3==1)
test5=subset(test4,B4==0 | B4==1)

founders=c(15:18)
test5$sum=apply(test5[,founders],1,sum)

test6=subset(test5,test5$sum != 0)
test7=subset(test6,test6$sum != 4)

index=c(1:4,15:19)
test8=test7[index]
founder_states=test8

write.table(founder_states,file="founder_states.txt",quote=FALSE)

