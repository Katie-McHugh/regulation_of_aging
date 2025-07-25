#############################################################################
### Budscars Plot
#############################################################################

## Load Libraries
#------------------------------------------------------------------------------
library(writexl)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(DescTools)
library(tidyverse)
library(magrittr)
library(rstatix)
#------------------------------------------------------------------------------

#read in and manage data
#------------------------------------------------------------------------------
scars2=read.table("data/clean/scarsdata.txt",header=T)
scars2$strain=as.factor(scars2$strain)
scars=subset(scars2, replicate!="2.7")
scars$strain <- as.character(scars$strain)
scars$strain[scars$strain == '4x'] <- '4S'
scars$strain <- as.factor(scars$strain)
#------------------------------------------------------------------------------

#############################################################################
### Statistics
#############################################################################

## Table 1
#------------------------------------------------------------------------------
### Calculating Averages for young and old fractions

scarsold=subset(scars, age=="aged")
scarsyoung=subset(scars, age=="young")

SO4S=subset(scarsold, strain=="4S")
SOyps128=subset(scarsold, strain=="YPS128")
SOy12=subset(scarsold, strain== "Y12")
SOdbvpg6044=subset(scarsold, strain=="DBVPG6044")
SOdbvpg6765=subset(scarsold, strain=="DBVPG6765")

SY4S=subset(scarsyoung, strain=="4S")
SYyps128=subset(scarsyoung, strain=="YPS128")
SYy12=subset(scarsyoung, strain=="Y12")
SYdbvpg6044=subset(scarsyoung, strain=="DBVPG6044")
SYdbvpg6765=subset(scarsyoung, strain=="DBVPG6765")

#old
summary(SO4S$avg_scars_all)
summary(SOyps128$avg_scars_all)
summary(SOy12$avg_scars_all)
summary(SOdbvpg6044$avg_scars_all)
summary(SOdbvpg6765$avg_scars_all)

#young
summary(SY4S$avg_scars_all)
summary(SYyps128$avg_scars_all)
summary(SYy12$avg_scars_all)
summary(SYdbvpg6044$avg_scars_all)
summary(SYdbvpg6765$avg_scars_all)

#------------------------------------------------------------------------------
### Checking normality of residuals
#------------------------------------------------------------------------------

scars4S=subset(scars, strain=="4S")
scarsY12=subset(scars, strain=="Y12")
scarsYPS128=subset(scars, strain=="YPS128")
scars6044=subset(scars, strain=="DBVPG6044")
scars6765=subset(scars, strain=="DBVPG6765")

scars4S_lm=lm(scars4S$avg_scars_all~scars4S$age)
scars4S_lm
#residuals
res.scars4S<-rstandard(scars4S_lm)
res.scars4S
shapiro.test(res.scars4S) #p=0.7092 #Not significant #t-test

scarsYPS128_lm=lm(scarsYPS128$avg_scars_all~scarsYPS128$age)
scarsYPS128_lm
#residuals
res.scarsYPS128<-rstandard(scarsYPS128_lm)
res.scarsYPS128
shapiro.test(res.scarsYPS128) #p=0.976 #Not significant #t-test

scarsY12_lm=lm(scarsY12$avg_scars_all~scarsY12$age)
scarsY12_lm
#residuals
res.scarsY12<-rstandard(scarsY12_lm)
res.scarsY12
shapiro.test(res.scarsY12) #p=0.00284 #significant #WSR

scars6044_lm=lm(scars6044$avg_scars_all~scars6044$age)
scars6044_lm
#residuals
res.scars6044<-rstandard(scars6044_lm)
res.scars6044
shapiro.test(res.scars6044) #p=0.173 #not significant #t-test

scars6765_lm=lm(scars6765$avg_scars_all~scars6765$age)
scars6765_lm
#residuals
res.scars6765<-rstandard(scars6765_lm)
res.scars6765
shapiro.test(res.scars6765) #p=0.5655 #not significant #t-test

#------------------------------------------------------------------------------
### Performing Pairwise Statistical Comparisons
#------------------------------------------------------------------------------

compare_means(avg_scars_all~age, data=scars, group.by="strain", paired=TRUE, method="t.test")
# 4x: p=0.00389
# YPS128: 0.00091
# Y12: Use WSR
# DBVPG6044: 0.00453
# DBVPG6765: 0.0284

compare_means(avg_scars_all~age, data=scars, group.by="strain", paired=TRUE)
# y12: p=0.0625

#------------------------------------------------------------------------------

#############################################################################
### Plotting
#############################################################################

## Figure 1
#------------------------------------------------------------------------------

#read in and manage data
# scars2=read.table("scarsdata.txt",header=T)
# scars2$strain=as.factor(scars2$strain)
# scars=subset(scars2, replicate!="2.7")
# scars$strain <- as.character(scars$strain)
# scars$strain[scars$strain == '4x'] <- '4S'
# scars$strain <- as.factor(scars$strain)
# View(scars)

library(viridis)
#set Theme
My_Theme = theme(
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16), 
  axis.text.y = element_text(size = 14))

# Use the Okabe-Ito colorblind-friendly palette
okabe_ito <- palette.colors(8, palette = "Okabe-Ito")
custom_colors <- c("young" = okabe_ito[7], "aged" = okabe_ito[6])

#bud scars plot as PDF
pdf("figures/SuppFig1_BUDSCARS_v2.pdf",height=4, width=8)
scars$strain<-factor(scars$strain, c("4S", "YPS128" ,"Y12", "DBVPG6044", "DBVPG6765"))
scarsplot<-ggplot(scars, aes(x = age, y = avg_scars_all, group = replicate))+
  geom_point(aes(color = age), size = 3)+
  geom_line()+
  labs(y = "Average bud scars per cell")+
  facet_wrap(~ strain, 
             strip.position = "top", ncol = 5)+   #put facet label to the top 
  theme_pubr()+
  scale_color_manual(values = custom_colors) +
  #scale_color_manual(values = c("young" = "orange", "aged" = "steelblue4"))+
  theme(strip.placement = "outside",   #move the facet label
        strip.text = element_text(size = 14),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1))
scarsplot+scale_x_discrete(limits=c("young", "aged")) +My_Theme
dev.off()

