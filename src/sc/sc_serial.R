# ==============================================================================
# Single Cell Samples Analysis
#
# Author: Max & Klara
#
# Description: Seq Analysis of SC samples and follow up
#
# Input: seqdata
#
# Output: plots
#
# ==============================================================================
########   Dependencies   #####
library(base)
library(dplyr)
library(xlsx)
library(stringr)
library(reshape)
library(tidyr)
library(readxl)
library(reshape2)
library(ggplot2)
library(ggsci)


########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')

load("data/interim/newsamples_filtered.RData")
load('data/interim/seqdata.RData')
SC_registry <- read.xlsx("data/external/SC_registry.xlsx",sheetIndex=1,header=TRUE)


################### Overview table ####
#works once tags are included in seqdata preparation
left_join(df,SC_registry, by="Sample_orig")%>%
  mutate(Patient.ID = Patient.ID.y)%>%
  mutate(Visite=Visite.y)%>%
  filter(Patient.ID %in% SC_registry$Patient.ID) ->df.sc

df.sc %>% 
  dplyr::select(Patient.ID,Sample_orig) %>% 
  unique %>% 
  dplyr::select(Patient.ID) %>% 
  table %>% 
  data.frame -> no_samples

df.sc %>% 
  left_join(.,no_samples,by="Patient.ID") %>%
  filter(Freq > 1)%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  group_by(Patmut) %>%
  mutate(n.mut = n()) %>%
  data.frame %>%
  filter(n.mut > 1) %>%
  filter(p.binom < -12) %>%
  filter(AF < 0.001) %>%
  filter(!snp)%>%
  filter(Func!="intronic")%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.01,
         minVAF < 0.38)-> serial.mut

df.sc %>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  filter(Patmut %in% serial.mut$Patmut)%>%
  filter(Sample_orig!="SC-3_FU")%>%  ###we don't want this in the plot
  mutate(timepoint = ifelse(Visite == 1, 0, tFU))%>%
ggplot() + 
  geom_point(aes(x=timepoint,y=TVAF,color=Gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=TVAF,group=position,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=3, dir="h") +
  scale_y_log10() +
  scale_color_igv()+
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.serial

ch_genes <- c("DNMT3A","TET2","ASXL1","PPM1D","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5")
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")

##Patient.ID
select(df.filtered,Patient.ID)%>%
  unique()->Patient.ID



