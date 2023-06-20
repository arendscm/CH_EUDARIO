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
library(viridis)
library(ggpubr)
library(g3viz)
library(ggsci)

#####  Data preparation ####
########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/clin.RData')
df -> df.clin
load('data/interim/seqdata.RData')
load('data/interim/seqdata_filtered.RData')
load('data/interim/seqdata_filtered_cf.RData')


######## Global input
source("src/material_table.R") #patient ids
source("src/genegroup_definitions.R")
######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")

##samples failed in sequencing
failedSamples <-c('OvCA_44_C1D1_cf','OvCA_45_C1D1_cf','OvCA_46_C1D1_cf','OvCA_48_C1D1_cf','OvCA_50_C1D1_cf','OvCA_54_C1D1_cf','OvCA_93_C1D1_cf',
                  'OvCA_11_C1D1_cf','OvCA_40_C1D1_cf','OvCA_53_C1D1_cf','OvCA_65_C1D1_cf')

##relevant variables
variables <- c("Patient.ID","Sample_orig","mutID","position","Sample", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "GeneDetail", "ExonicFunc", "AAChange", "cytoBand","readDepth", "TR1", "TR1_plus", "TR1_minus", "TR2", "TR2_plus", "TR2_minus", "TVAF", "AF", "avsnp150","cosmic92_coding","snp","mutFreq","p.binom","n.mut","n.material","sum_cf","sum_wb","Material","tag", "Patmut")


#####  serial analysis: cf data exploration 2 timepoints cf only ####
df %>% 
  filter(Material=="cf")%>%
  filter(c1d1_cf==1&eot_cf==1)%>%  
  #filter(c1d1_cf==1&eot_cf==1 & c7d1_cf == 0)%>% #for subsetting patients with 2 or 3 timepoints
  filter(n.material>=2)%>% 
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.01)%>%
  filter(snp==FALSE)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes_without_HRD),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA","other")))))%>%
  filter(maxVAF > 0.005,
         maxTR2 > 7,
         minpbinom < -8) %>%
  #filter(TVAF < 0.38) %>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=position,color=gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=8, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.cf.serial

png("output/figures/p.cf.serial.png",width=10, height=5,units="in",res=500,type="cairo")
p.cf.serial
dev.off()

#####  serial analysis: data preparation: find mutations that are present in both wb and cf  #####
df%>%  mutate(cfID=paste(Patient.ID,position,sep="_"))%>%
  filter(Material=="wb")%>%
  group_by(Visite)%>%
  dplyr::select(cfID)->wb.cfID

df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(Material=="cf")%>%
  group_by(Visite)%>%
  mutate(wb=is.element(mutID,wb.mutID))->cf.cfID

#####  serial analysis: cf data exploration 2 timepoints cf and wb - wb=full line ; plasma=dashed ####
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=gene,group=pat_material,shape=Material),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,color=gene,linetype=Material),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=6, scales="free", dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial


#####  Serial analysis: facet by patient and material including "cf_only" ####
##serial cf data exploration 2 timepoints cf and wb
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n()) %>% data.frame()%>% filter(n.visit==1)%>% filter(Material=="cf") %>% mutate(Material="cf_only") ->df.cfonly

df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  full_join(.,df.cfonly)%>%
  filter(is.element(Patient.ID,c("4202011","4207004","4202001","4221010")))%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  filter(Patient.ID %in% c("4203009","4203002", "4221004", "4221006", "4221023"))%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=gene,group=pat_material),size=1,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,color=gene),size=0.5,na.rm=FALSE) + 
  facet_grid(Material ~ Patient.ID) +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial

png("output/figures/p.cf.serial_cf,cf_only,wb_Patientsexample.png",width=6, height=4,units="in",res=500,type="cairo")
p.cf.serial
dev.off()


#####  Serial analysis: single patient by mutation ####
df %>% 
  filter(Patient.ID=="4221010")%>%
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes_without_HRD),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  mutate(AAChange_mod = str_replace(AAChange, ":NM.*:p.", ":p."))%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=Gene,group=pat_material,shape=Material),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,color=Gene,linetype=Material),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ AAChange_mod, ncol=6, dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  #scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial

##serial cf data exploration 2 timepoints cf and wb TP53 only
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(Gene=="TP53")%>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding,"ovary")) %>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  filter(maxVAF > 0.001,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=factor(n.visit),group=pat_material,shape=Material),size=2,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,linetype=Material),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=5, dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="detected in WB/cfDNA") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial.tp53


##### serial cf data exploration 2 timepoints cf and wb HRD only
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(is.element(Gene,hrd_genes))%>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding,"ovary")) %>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  filter(maxVAF > 0.005,
         maxTR2 > 19,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=factor(n.visit),group=pat_material,shape=Material),size=2,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,linetype=Material,color=cosmic_ovary),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=5, dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="detected in WB/cfDNA") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial


###test cf serial with BRCA mutations
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding,"ovary")) %>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  filter(gene=="BRCA")%>%
  filter(maxVAF > 0.005,
         maxTR2 > 19,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=factor(n.visit),group=pat_material,shape=Material),size=2,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,linetype=Material,color=cosmic_ovary),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=5, dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="detected in WB/cfDNA") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()

#####  just for fun: dynamics cfVAF vs wbVAF ####
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  filter(Material == "wb")%>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))->df1
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  filter(Material == "cf")%>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))->df2

full_join(df1,df2,by=c("cfID_visit"))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  filter(TVAF.x!=0)%>%
  ggplot(aes(x=log(TVAF.x+0.001),y=log10(TVAF.y+0.001),
             color=gene.x)) +
  geom_point(aes(shape=Visite.y),size=3)+
  geom_line(aes(group=Patmut.y),size=1)+
  #scale_y_log10()+
  #scale_x_log10()+
  scale_color_viridis(discrete=TRUE)+
  ylab("whole blood VAF")+
  xlab("cfDNA VAF")+
  theme_minimal()

#png("output/figures/p.cf.serial.png",width=6, height=4,units="in",res=500,type="cairo")
#p.cf.serial
#dev.off()


df %>% 
  filter(Material=="cf")%>%
  filter(c1d1_cf==1&c7d1_cf==1&eot_cf==1)%>%
  filter(Patient.ID %in% c("4203001", "4203004"))%>%
  filter(n.material==3)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.01)%>%
  filter(snp==FALSE)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA","other")))))%>%
  filter(maxVAF > 0.005,
         maxTR2 > 7,
         minpbinom < -8)%>%
  .$Patmut->Patmut.serial

df%>%
  filter(Patmut %in% Patmut.serial)%>%
  filter(Material == "cf")%>%
  #filter(TVAF < 0.38) %>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=Gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=position,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=2, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.3)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal() -> p.cf.serial

png("output/figures/p.cf.serialC7.png",width=6, height=4,units="in",res=500,type="cairo")
p.cf.serial
dev.off()


####### serial analysis in patients with 3 timepoints 
df %>% 
  filter(!is.na(Patient.ID))%>%
  filter(is.na(replicate))%>%
  filter(c1d1_wb==1) %>% 
  filter(c1d1_cf==1&c7d1_cf==1&(ue_cf==1|eot_cf==1))%>%
  filter(sum_cf==3)%>%
  filter(Material=="cf")%>%
  mutate(timepoint = ifelse(Visite == "C1D1",0,
                            ifelse(Visite == "EOT",tEOT,
                                   ifelse(Visite=="UE",tUE,
                                          ifelse(Visite=="C7D1",tC7D1,NA)))))-> df.cf_serial

df.cf_serial %>% 
  filter(n.material>2)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(p.binom < -10) %>% 
  mutate(fraction_mut = mutFreq/n.lane) %>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.008) %>%
  filter(minVAF < 0.35) %>%
  filter(fraction_mut < 0.2) %>%
  .$Patmut -> Patmut.serial

###all clones by patient
df.cf_serial%>%
  filter(Patmut %in% Patmut.serial)%>%
  filter(Material == "cf")%>%
  #filter(Gene %in% ch_genes)%>%
  #filter(TVAF < 0.38) %>%
  ggplot() + 
  geom_point(aes(x=timepoint,y=TVAF,color=Gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=TVAF,group=Patmut,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=6, dir="h") +
  scale_y_log10()+
  scale_color_igv()+
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal() -> p.cf.serial


df.cf_serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "C1D1") %>%
  mutate(vaf_c1d1=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_c1d1,timepoint)-> df.cf_c1d1

df.cf_serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "C7D1") %>%
  mutate(vaf_c7d1=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_c7d1,timepoint)-> df.cf_c7d1

df.cf_serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "UE") %>%
  mutate(vaf_eot=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_eot,timepoint)-> df.cf_ue

df.cf_serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "EOT") %>%
  mutate(vaf_eot=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_eot,timepoint)-> df.cf_eot

full_join(df.cf_ue,df.cf_eot) -> df.cf_eotue

df.cf_c1c7 <- full_join(df.cf_c1d1,df.cf_c7d1,by=c("Patient.ID","Gene","AAChange","ExonicFunc","position")) %>% 
  mutate(timepoint = timepoint.y-timepoint.x,
         relvaf1 = vaf_c1d1/vaf_c1d1,
         relvaf2 = vaf_c7d1/vaf_c1d1) %>%
  mutate(fitness = log(vaf_c7d1/vaf_c1d1*(1-2*vaf_c1d1)/(1-2*vaf_c7d1))/(timepoint/365))%>%
  melt.data.frame(measure.vars = c("relvaf1","relvaf2"))

df.cf_c7eot <- full_join(df.cf_c7d1,df.cf_eotue,by=c("Patient.ID","Gene","AAChange","ExonicFunc","position")) %>% 
  mutate(timepoint = timepoint.y-timepoint.x,
         relvaf1 = vaf_c7d1/vaf_c7d1,
         relvaf2 = vaf_eot/vaf_c7d1) %>%
  mutate(fitness = log(vaf_eot/vaf_c7d1*(1-2*vaf_c7d1)/(1-2*vaf_eot))/(timepoint/365))%>%
  melt.data.frame(measure.vars = c("relvaf1","relvaf2"))


###direct comparison of clone fitnesss c1c7 vs c7eot only in patients which had 5 or 6 cycles of Platinum 
full_join(df.cf_c1c7 %>% filter(variable=="relvaf2")%>%select(Patient.ID,position,Gene,fitness),
          df.cf_c7eot%>%filter(variable=="relvaf2")%>%select(Patient.ID,position,Gene,fitness),
          by=c("Patient.ID","Gene","position")) %>%
  filter(Patient.ID %in% (df.clin %>% filter(Number_ChemotherapyCycles > 5))$Patient.ID)%>%
  mutate(fitness_platinum = fitness.x, 
         fitness_maintenance = fitness.y)%>%
  dplyr::select(-fitness.x,-fitness.y)%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  melt.data.frame(measure.vars = c("fitness_platinum","fitness_maintenance"))%>%
  filter(Gene %in% c("CHEK2","PPM1D","DNMT3A","TP53","TET2"))%>%
  ggplot(aes(x=variable,y=value,color=Gene))+
  geom_boxplot(color="grey")+
  geom_point()+
  geom_line(aes(group=Patmut))+
  ylim(c(-15,15))+
  scale_color_npg()+
  facet_grid(~Gene)+
  theme_minimal()

  
#### boxplot with p values

full_join(df.cf_c1c7 %>% filter(variable=="relvaf2")%>%select(Patient.ID,position,Gene,fitness),
          df.cf_c7eot%>%filter(variable=="relvaf2")%>%select(Patient.ID,position,Gene,fitness),
          by=c("Patient.ID","Gene","position")) %>%
  filter(Patient.ID %in% (df.clin %>% filter(Number_ChemotherapyCycles > 5))$Patient.ID)%>%
  mutate(fitness_platinum = fitness.x, 
         fitness_maintenance = fitness.y)%>%
  dplyr::select(-fitness.x,-fitness.y)%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  melt.data.frame(measure.vars = c("fitness_platinum","fitness_maintenance"))%>%
  filter(Gene %in% c("CHEK2","PPM1D","DNMT3A","TP53","TET2"))%>%
  ggboxplot(.,
              x = "variable",
              y = "value",
              #combine = TRUE,
              #color = "DDR", 
              # palette = ,
              facet.by= "Gene",
              xlab = "Therapy",
              ylab = "Fitness",
              title = "",
              width = 0.3,
              #ylim = c(-5,10),
              size=0.8,
              alpha=1,
              repel=TRUE,
              #yscale = "log10",
              scales = "free",
              add = c("jitter")
  )+  stat_compare_means()


##scatterplot growthrate Platinum+PARPi vs PARPI maintenance
full_join(df.cf_c1c7 %>% filter(variable=="relvaf2")%>%select(Patient.ID,position,Gene,fitness),
          df.cf_c7eot%>%filter(variable=="relvaf2")%>%select(Patient.ID,position,Gene,fitness),
          by=c("Patient.ID","Gene","position")) %>%
  filter(Patient.ID %in% (df.clin %>% filter(Number_ChemotherapyCycles > 5))$Patient.ID)%>%
  mutate(fitness_platinum = fitness.x, 
         fitness_maintenance = fitness.y)%>%
  dplyr::select(-fitness.x,-fitness.y)%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  group_by(Gene)%>%
  summarise(med.fit.plat = mean(fitness_platinum),
            med.fit.maint = mean(fitness_maintenance),
            err.fit.plat = sd(fitness_platinum)/sqrt(n()),
            err.fit.maint = sd(fitness_maintenance)/sqrt(n()),
            size=n())%>%
  data.frame %>%
  filter(size>5)%>%
  ggplot(aes(x=med.fit.plat,y=med.fit.maint,color=Gene))+
  geom_point(aes(size=size),alpha = 1)+
  geom_errorbar(aes(xmin=med.fit.plat-err.fit.plat,xmax=med.fit.plat+err.fit.plat),width=0.2)+
  geom_errorbar(aes(ymin=med.fit.maint-err.fit.maint,ymax=med.fit.maint+err.fit.maint),width=0.2)+
  scale_color_npg()+
  geom_abline(slope=1,intercept=0)+
  ylim(c(-5,5))+
  xlim(c(-5,5))
  theme_minimal()
  
  
  
  
  
  ######### identity check for sequential plasma samples
  df %>% 
    filter(!is.na(Patient.ID))%>%
    filter(is.na(replicate))%>%
    filter(c1d1_wb==1) %>% 
    filter(c1d1_cf==1&c7d1_cf==1&(ue_cf==1|eot_cf==1))%>%
    filter(sum_cf==3)%>%
    filter(Material=="cf")%>%
    mutate(timepoint = ifelse(Visite == "C1D1",0,
                              ifelse(Visite == "EOT",tEOT,
                                     ifelse(Visite=="UE",tUE,
                                            ifelse(Visite=="C7D1",tC7D1,NA)))))-> df.cf_serial
  
  
  df.cf_serial %>% 
    filter(snp)%>%
    filter(Visite == "C1D1") %>%
    mutate(vaf_c1d1=TVAF) %>%
    dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_c1d1,timepoint)-> df.cf_c1d1
  
  df.cf_serial %>% 
    filter(snp)%>%
    filter(Visite == "C7D1") %>%
    mutate(vaf_c7d1=TVAF) %>%
    dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_c7d1,timepoint)-> df.cf_c7d1
  
  df.cf_serial %>% 
    filter(snp)%>%
    filter(Visite == "UE") %>%
    mutate(vaf_eot=TVAF) %>%
    dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_eot,timepoint)-> df.cf_ue
  
  df.cf_serial %>% 
    filter(snp)%>%
    filter(Visite == "EOT") %>%
    mutate(vaf_eot=TVAF) %>%
    dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_eot,timepoint)-> df.cf_eot
  
  
  join_var = c("Patient.ID","Gene","AAChange","ExonicFunc","position")
  full_join(df.cf_c1d1,
            full_join(df.cf_c7d1,
                      full_join(df.cf_ue,df.cf_eot,by=join_var),by=join_var),by=join_var)