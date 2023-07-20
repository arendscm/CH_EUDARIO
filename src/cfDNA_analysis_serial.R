# ______________________________________________________________________________
# Longitudinal analysis in cfDNA and WB
#
# Author: Max & Klara
#
# Description: Analysis of longitudinal data cfDNA samples
#
# Input: seqdata
#
# Output: serial sample analysis of cf DNA and wb DNA mutations
#
# ______________________________________________________________________________
#####  Dependencies   #####
library(base)
library(dplyr)
library(stringr)
library(reshape)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(g3viz)
library(ggsci)
library(survminer)
library(survival)

#####  Data preparation ####
########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/clin.RData')
load('data/interim/seqdata.RData')
load('data/interim/seqdata_filtered.RData')
load('data/interim/seqdata_filtered_cf.RData')


######## Global input
source("src/material_table.R") #patient ids
source("src/genegroup_definitions.R")
######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")

##relevant variables
variables <- c("Patient.ID","Sample_orig","mutID","position","Sample", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "GeneDetail", "ExonicFunc", "AAChange", "cytoBand","readDepth", "TR1", "TR1_plus", "TR1_minus", "TR2", "TR2_plus", "TR2_minus", "TVAF", "AF", "avsnp150","cosmic92_coding","snp","mutFreq","p.binom","n.mut","n.material","sum_cf","sum_wb","Material","tag", "Patmut")


####### serial analysis in patients with 2 timepoints c1d1 and c7d1
df %>% 
  filter(is.na(replicate))%>%
  filter(!is.na(Patient.ID))%>%
  filter(c1d1_cf==1&c7d1_cf==1)%>%
  filter(Material=="cf")%>%
  filter(Visite=="C1D1"|Visite=="C7D1")%>%
  mutate(timepoint = ifelse(Visite == "C1D1",0,
                            ifelse(Visite=="C7D1",tC7D1,NA)))-> df.cf_serial


df.cf_serial %>% 
  group_by(Patmut) %>%
  mutate(n.serial = n())%>%
  data.frame%>%
  filter(n.serial > 1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.01)%>%
  filter(snp==FALSE)%>%
  filter(p.binom < -12) %>% 
  mutate(fraction_mut = mutFreq/n.lane) %>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         max_n.visite = max(n.visite)) %>%
  data.frame()%>%
  filter(max_n.visite == 1)%>% # this should restrict to cf only clones
  filter(maxVAF > 0.005) %>%
  filter(minVAF < 0.35) %>%
  filter(fraction_mut < 0.2) %>%
  filter(tag != "true")%>%
  #filter(Patient.ID %in% patients.bc)%>%
  .$Patmut -> Patmut.serial2

###include all mutations tagged true (some are sub threshold at certain timepoints)
df.cf_serial %>% 
  filter(tag=="cf-only")%>%  
  group_by(Patmut) %>%
  mutate(n.serial = n())%>%
  data.frame%>%
  .$Patmut -> Patmut.serial1

Patmut.serial <- c(Patmut.serial1,Patmut.serial2) %>% unique

####
df.cf_serial %>%
  filter(Patmut %in% Patmut.serial) %>%
  group_by(Patmut) %>%
  mutate(n.serial = n())%>%
  data.frame%>%
  filter(n.serial == 1) %>%
  mutate(Visite = ifelse(Visite=="C1D1","C7D1","C1D1"),
         TVAF = 0.001,
         timepoint = ifelse(Visite=="C1D1",0,tC7D1)) -> df.missing  ##manually changes the visit and timepoint and sets VAF to 0.001 (threshold) at the missing timepoint

df.missing$Patmut -> Patmut_missing

df.cf.serial <- full_join(df.cf_serial,df.missing)
###all clones by patient

df.cf.serial%>%
  filter(Patmut %in% Patmut.serial)%>%
  filter(Gene %in% c("TP53","NF1","CDK12","EMSY"))%>%
  ggplot() + 
  geom_point(aes(x=timepoint,y=TVAF,color=Gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=TVAF,group=Patmut,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=10, dir="h") +
  scale_y_log10()+
  scale_color_igv()+
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal() -> p.cf.serial


df.cf.serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "C1D1") %>%
  mutate(vaf_c1d1=TVAF,
         depth_c1=readDepth) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_c1d1,timepoint,depth_c1)-> df.cf_c1d1

df.cf.serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "C7D1") %>%
  mutate(vaf_c7d1=TVAF,
         depth_c7=readDepth) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_c7d1,timepoint,depth_c7)-> df.cf_c7d1


df.cf_c1c7 <- full_join(df.cf_c1d1,df.cf_c7d1,by=c("Patient.ID","Gene","AAChange","ExonicFunc","position")) %>% 
  mutate(timepoint = timepoint.y-timepoint.x,
         sig.x = sqrt((vaf_c1d1*(1-vaf_c1d1))/depth_c1), 
         sig.y = sqrt((vaf_c7d1*(1-vaf_c7d1))/depth_c7),
         xy = ifelse((vaf_c7d1 < vaf_c1d1-1.96*sig.x)|(vaf_c7d1 > vaf_c1d1+1.96*sig.x),0,1),
         yx = ifelse((vaf_c1d1 < vaf_c7d1-1.96*sig.y)|(vaf_c1d1 > vaf_c7d1+1.96*sig.y),0,1),
         relvaf1 = vaf_c1d1/vaf_c1d1,
         relvaf2 = vaf_c7d1/vaf_c1d1) %>%
  
  melt.data.frame(measure.vars = c("relvaf1","relvaf2"))



###direct comparison of clone fitnesss c1c7 vs c7eot only in patients which had 5 or 6 cycles of Platinum 
df.cf_c1c7 %>%
  filter(Patient.ID %in% (df.clin %>% filter(Number_ChemotherapyCycles > 5))$Patient.ID)%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  mutate(timepoint = ifelse(variable=="relvaf1",0,timepoint))%>%
  mutate(fitness_binom = ifelse(xy==1|yx==1,"stable",
                                ifelse(vaf_c7d1/vaf_c1d1 < 1,"decreasing","increasing")))%>%
  mutate(category = ifelse(xy==1|yx==1,0,
                           ifelse(vaf_c7d1/vaf_c1d1<0.1,-1,
                           ifelse(vaf_c7d1/vaf_c1d1>10,1,0))))%>%
  filter(Gene == "TP53")%>%
  filter(timepoint.y>50)%>%
  ggplot(aes(x = timepoint, y = value, color = factor(category), group = Patmut)) +
  geom_point(size = 1.5, na.rm = FALSE) + 
  geom_line()+
  scale_y_log10() +
  scale_color_manual(values=c("#4DBBD5FF","#E64B35FF","#00A087FF"))+
  labs(x = "Time in days", y = "VAF change", color = "Fitness") +
  my_theme() +
  theme(strip.text = element_text(face = "italic"))-> p.growth
p.growth

###correlation with response

df.cf_c1c7 %>%
  filter(Patient.ID %in% (df.clin %>% filter(Number_ChemotherapyCycles > 5))$Patient.ID)%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  mutate(timepoint = ifelse(variable=="relvaf1",0,timepoint))%>%
  mutate(fitness_binom = ifelse(xy==1|yx==1,0,
                                ifelse(vaf_c7d1/vaf_c1d1 < 1,-1,1)))%>%
  mutate(category = ifelse(xy==1|yx==1,0,
                           ifelse(vaf_c7d1/vaf_c1d1<0.1,-1,
                                  ifelse(vaf_c7d1/vaf_c1d1>10,1,0))))%>%
  filter(timepoint.y>50)%>%
  filter(Gene %in% c("TP53","NF1","CDK12","EMSY"))%>%
  group_by(Patient.ID) %>%
  mutate(category_max = max(category))%>%
  data.frame %>% 
  dplyr::select(Patient.ID,category_max) %>% 
  unique %>%
  left_join(.,df.clin,by=c("Patient.ID")) -> df.clin_tp53


##assessing best response
df.clin_tp53 %>% dplyr::select(category_max,Response_best)%>% table 

df.clin_tp53 %>% 
  dplyr::select(category_max,Response_best)%>% 
  table%>%
  data.frame%>%
  ggplot(aes(x=category_max,y=Freq,fill=Response_best))+geom_bar(stat="identity")
##assessing PFS
df.clin_tp53 -> df.surv
surv_obj <- Surv(time = df.surv$PFS_days, event = df.surv$PFS_event)
fit.km <- survfit(surv_obj ~ category_max, data = df.surv)

df.surv %>% 
  ggsurvplot(fit.km, 
             data = ., 
             pval = TRUE,
             #conf.int = TRUE,
             palette = pal_npg("nrc")(3),
             risk.table=TRUE,
             tables.height = 0.3,
             ylab = "Progression-free Survival",
             xlab = "Time in Days",
             #legend.title = "CH status",
             #legend.labs = c("negative", "positive")
  ) -> p.pfs

p.pfs



surv_obj <- Surv(time = df.surv$OS_days, event = df.surv$OS_event)
fit.km <- survfit(surv_obj ~ category_max, data = df.surv)

df.surv %>% 
  ggsurvplot(fit.km, 
             data = ., 
             pval = TRUE,
             #conf.int = TRUE,
             palette = pal_npg("nrc")(3),
             risk.table=TRUE,
             tables.height = 0.3,
             ylab = "Overall Survival",
             xlab = "Time in Days",
             #legend.title = "CH status",
             #legend.labs = c("negative", "positive")
  ) -> p.os

p.os


##Correlation of TP53 VAF with CA125 and tumor volume
library(GGally)

df.cf_c1c7 %>% 
  filter(is.element(Gene,c("TP53")))%>%
  group_by(Patient.ID)%>%
  mutate(maxVAF = max(vaf_c1d1))%>%
  data.frame%>%
  filter(vaf_c1d1==maxVAF)%>%
  left_join(.,df.clin %>% dplyr::select(Patient.ID,CA125,TumorBurden_baseline),by="Patient.ID") %>% 
  mutate(logVAF = log10(vaf_c1d1),
         logCA = log10(CA125))%>%
  dplyr::select(logCA,TumorBurden_baseline,logVAF) %>% 
  ggpairs(.)
  
