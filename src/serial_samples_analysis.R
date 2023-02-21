# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: Analysis of longitudinal data
#
# Input: seqdata, seqdata_filtered
#
# Output: ...
#
# ==============================================================================

########   Dependencies   #####
library(base)
library(dplyr)
library(stringr)
library(reshape)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(viridis)
library(reshape)
library(ggpubr)
library(maftools)

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')

######## Get Patient ids
source("src/ids.R")

######## Functions and themes

source("src/global_functions_themes.R")

########   SERIAL SAMPLES ####

####  identity check via SNPs ####
df %>% 
  filter(Visite=="EOT") %>% 
  mutate(ID=paste(Patient.ID,position))-> df.eot1

ids %>% filter(Visite=="EOT") -> eotsamples

df %>% 
  filter(Visite=="C1D1") %>% 
  filter(is.element(Patient.ID,eotsamples$Patient.ID)) %>%
  mutate(ID=paste(Patient.ID,position)) -> df.c1d0

left_join(df.eot1,df.c1d0,by="ID") %>% 
  filter(snp.x == 1) %>% 
  ggplot(aes(x=factor(Patient.ID.x),y=TVAF.x-TVAF.y)) +
  geom_point()

##serial samples dynamics plot
df%>% 
  filter(!is.na(Patient.ID))%>%
  filter(is.na(replicate))%>%
  filter(c1d1_wb==1 & eot_wb==1)%>% 
  filter(Material=="wb")%>%
  filter(Visite == "C1D1")-> df.eot.c1
df%>% 
  filter(!is.na(Patient.ID))%>%
  filter(is.na(replicate))%>%
  filter(c1d1_wb==1 & eot_wb==1)%>% 
  filter(Material=="wb")%>%
  filter(Visite == "EOT")-> df.eot.eot
semi_join(df.eot.c1,df.eot.eot, by = "Patmut")->df.eot
df.eot$Patmut%>%
  unique()->Patmut.serial

df %>% 
  filter(!is.na(Patient.ID))%>%
  filter(is.na(replicate))%>%
  filter(c1d1_wb==1&eot_wb==1) %>% 
  filter(Material=="wb")%>%
  filter(Patmut %in% Patmut.serial)-> df.eot

df.eot %>% 
  filter(n.material>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(p.binom < -10) %>% 
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.008) %>%
  filter(TVAF < 0.38) %>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=position,color=Gene),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=6, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.serial

png("output/figures/p.serial.png",width=6, height=4,units="in",res=500,type="cairo")
p.serial
dev.off()

##Examples
df.eot %>% 
  filter(is.element(Patient.ID,c("4221010","4202008","4202016")))%>%
  filter(n.material>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(p.binom < -10) %>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.008) %>%
  filter(TVAF < 0.38) %>%
  filter(Gene!="CEBPA")%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=position,color=Gene),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=6, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Timepoint",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.serial

png("output/figures/p.serial.png",width=6, height=3,units="in",res=500,type="cairo")
p.serial
dev.off()
###Serial samples by brca status (question: do dynamics unter PARP Inhb. differ depending on BRCA status?)
source("src/brca_germline.R")
df.eot %>% 
  left_join(.,id.brca_germline,by = "Patient.ID")%>%
  filter(Gene=="PPM1D"|Gene=="TP53")%>%
  filter(serial.mut>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(mutFreq < 10) %>% 
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.008) %>%
  filter(TVAF < 0.38) %>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=position,color=Gene),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ brca_germline, ncol=2, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.serial

##TEST SERIAL SAMPLES relative >- for this we need the timedifference between d1 and eot

df.eot %>% 
  filter(serial.mut>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(mutFreq < 10) %>% 
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.005) %>%
  filter(TVAF < 0.37) %>%
  filter(Visite == "C1D1") %>%
  mutate(vaf_d1=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,position,vaf_d1)-> df.eotd1

df.eot %>% 
  filter(serial.mut>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(mutFreq < 10) %>% 
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.005) %>%
  filter(TVAF < 0.37) %>%
  filter(Visite == "EOT") %>%
  mutate(vaf_eot=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,position,vaf_eot)-> df.eoteot

df.eot_rel <- full_join(df.eotd1,df.eoteot) %>% mutate(relvaf1 = vaf_d1/vaf_d1,
                                                       relvaf2 = vaf_eot/vaf_d1) %>%
  melt.data.frame(measure.vars = c("relvaf1","relvaf2"))

df.eot_rel %>%  
  ggplot() + 
  geom_point(aes(x=variable,y=value,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=variable,y=value,group=position,color=Gene),size=1*1,na.rm=FALSE) + 
  scale_y_log10() +
  labs(x="Timepoint",y="log(VAF change)",colour="Mutated Gene") +
  theme_minimal()-> p.serial

df.eot_rel %>% 
  filter(variable == "relvaf2") %>% 
  ggplot() + 
  geom_point(aes(x=Gene,y=value,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  scale_y_log10() +
  labs(x="Gene",y="log(VAF change)",colour="Mutated Gene") +
  theme_graphicalabstract()-> p.serial
#hier stimmt was noch nicht, im plot ist eine lila linie vertikal bei relvaf2

df.eot_rel %>% 
  filter(variable == "relvaf2") %>% 
  filter(is.element(Gene,c("CHEK2","PPM1D","DNMT3A","TP53","TET2","ASXL1")))%>%
  mutate(DDR = ifelse(is.element(Gene,c("TP53","PPM1D","CHEK2")),"DDR","non-DDR"))%>%
  ggboxplot(., 
            x = "Gene",
            y = "value",
            #facet.by = "variable",
            # panel.labs = list(CHIP = c("positive","negative"), variable=c("Troponin","VCAM","hsCRP")), 
            combine = TRUE,
            color = "DDR", 
            palette = viridis(3),
            xlab = "Gene",
            ylab = "n-fold change in VAF",
            title = "",
            width = 0.3,
            ylim = c(0,15),
            size=0.8,
            alpha=1,
            repel=TRUE,
            #yscale = "log10",
            scales = "free",
            add = c("jitter"))+
  theme_minimal() + 
  theme(axis.title.x = element_blank()) +
  theme(#legend.position = "none",
        axis.text.x = element_text(face="italic"),
        axis.title.y = element_text(face ="plain"),
        plot.title = element_text(hjust=0,face ="plain")) ->p.growth

png("output/figures/growth.png",width=6, height=4,units="in",res=500,type="cairo")
p.growth
dev.off()


rm(list=ls())
