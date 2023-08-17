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

####### serial analysis in patients with 3 timepoints 
df %>% 
  filter(!is.na(Patient.ID))%>%
  filter(is.na(replicate))%>%
  filter(c1d1_wb==1) %>% 
  filter(c1d1_cf==1&c7d1_cf==1&(ue_cf==1|eot_cf==1))%>%
  filter(Material=="cf")%>%
  mutate(timepoint = ifelse(Visite == "C1D1",0,
                            ifelse(Visite == "EOT",tEOT,
                                   ifelse(Visite=="UE",tUE,
                                          ifelse(Visite=="C7D1",tC7D1,NA)))))%>%
  filter(Visite=="C1D1"|Visite=="C7D1"|Visite==lastTimepoint)-> df.cf_serial


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
  filter(max_n.visite ==2)%>% # this should restrict to CH clones
  filter(maxVAF > 0.005) %>%
  filter(minVAF < 0.35) %>%
  filter(fraction_mut < 0.2) %>%
  #filter(Patient.ID %in% patients.bc)%>%
  .$Patmut -> Patmut.serial2

###include all mutations tagged true (some are sub threshold at certain timepoints)
df.cf_serial %>% 
  filter(tag=="true")%>%  
  group_by(Patmut) %>%
  mutate(n.serial = n())%>%
  data.frame%>%
  filter(n.serial !=1)%>%
  .$Patmut -> Patmut.serial1

Patmut.serial <- c(Patmut.serial1,Patmut.serial2) %>% unique

####
df.cf_serial %>%
  filter(Patmut %in% Patmut.serial) %>%
  group_by(Patmut) %>%
  mutate(n.serial = n())%>%
  data.frame%>%
  filter(n.serial ==2) %>%
  mutate(tLTP = ifelse(lastTimepoint=="UE",tUE,tEOT))%>%
  mutate(missing.dummy = ifelse(Visite=="C1D1",1,
                                ifelse(Visite=="C7D1",2,4)))%>%
  group_by(Patmut)%>%
  mutate(missing = sum(missing.dummy),       ###this identifies the "missing" timepoint 3 = EOT/UE, 5 = C7D1, 6 = C1D1
         min.missing.dummy = min(missing.dummy))%>%   
  data.frame%>%
  filter(missing.dummy==min.missing.dummy)%>%  ##filters out one of the two time points
  mutate(Visite = ifelse(missing==3,"EOT",
                         ifelse(missing==5,"C7D1","C1D1")),
         TVAF = 0.001,
         timepoint = ifelse(Visite=="C1D1",0,
                            ifelse(Visite=="C7D1",tC7D1,tLTP))) %>%   ##manually changes the visit and timepoint and sets VAF to 0.001 (threshold) at the missing timepoint
  dplyr::select(-c(min.missing.dummy,missing.dummy,n.serial)) -> df.missing


df.cf.serial <- full_join(df.cf_serial,df.missing)
###all clones by patient

df.cf.serial%>%
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


df.cf.serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "C1D1") %>%
  mutate(vaf_c1d1=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_c1d1,timepoint)-> df.cf_c1d1

df.cf.serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "C7D1") %>%
  mutate(vaf_c7d1=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_c7d1,timepoint)-> df.cf_c7d1

df.cf.serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "UE") %>%
  mutate(vaf_eot=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_eot,timepoint)-> df.cf_ue

df.cf.serial %>% 
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
full_join(df.cf_c1c7 %>% filter(variable=="relvaf2")%>%dplyr::select(Patient.ID,position,Gene,fitness),
          df.cf_c7eot%>%filter(variable=="relvaf2")%>%dplyr::select(Patient.ID,position,Gene,fitness),
          by=c("Patient.ID","Gene","position")) %>%
  filter(Patient.ID %in% (df.clin %>% filter(Number_ChemotherapyCycles > 5))$Patient.ID)%>%
  mutate(fitness_platinum = fitness.x, 
         fitness_maintenance = fitness.y)%>%
  dplyr::select(-fitness.x,-fitness.y)%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  melt.data.frame(measure.vars = c("fitness_platinum","fitness_maintenance"))%>%
  filter(Gene %in% c("CHEK2","PPM1D","DNMT3A","TP53","TET2"))%>%
  ggplot(aes(x=variable,y=value,color=Gene))+
  geom_boxplot(color="darkgrey")+
  geom_point()+
  geom_line(aes(group=Patmut),alpha=0.5)+
  scale_x_discrete(labels=c("C","M"),name="Therapy")+
  scale_y_continuous(name="Fitness",limits=c(-10,15))+
  scale_color_npg()+
  facet_grid(~Gene)+
  my_theme() +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none",
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())-> p.fitness_CvsM
p.fitness_CvsM

png("output/figures/p.cf.serial_CvsM.png",width=6, height=3,units="in",res=500,type="cairo")
p.fitness_CvsM
dev.off()

#### boxplot with p values
my_comp =c("fitness_platinum","fitness_maintenance")
full_join(df.cf_c1c7 %>% filter(variable=="relvaf2")%>%dplyr::select(Patient.ID,position,Gene,fitness),
          df.cf_c7eot%>%filter(variable=="relvaf2")%>%dplyr::select(Patient.ID,position,Gene,fitness),
          by=c("Patient.ID","Gene","position")) %>%
  filter(Patient.ID %in% (df.clin %>% filter(Number_ChemotherapyCycles > 5))$Patient.ID)%>%
  mutate(fitness_platinum = fitness.x, 
         fitness_maintenance = fitness.y)%>%
  dplyr::select(-fitness.x,-fitness.y)%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  melt.data.frame(measure.vars = c("fitness_platinum","fitness_maintenance"))%>%
  filter(Gene %in% c("TP53","PPM1D"))%>%
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
              ylim = c(-10,15),
              size=0.8,
              alpha=1,
              repel=TRUE,
              #yscale = "log10",
              scales = "free",
              add = c("jitter")
  )+  
  scale_x_discrete(labels=c("C","M"),name="Therapy")+
  stat_compare_means(label.y =13,aes(label=paste0("p = ", after_stat(p.format))))+
  my_theme()+
  theme(strip.text = element_text(face = "italic"))
  


##scatterplot growthrate Platinum+PARPi vs PARPI maintenance mean+sem
full_join(df.cf_c1c7 %>% filter(variable=="relvaf2")%>%dplyr::select(Patient.ID,position,Gene,fitness),
          df.cf_c7eot%>%filter(variable=="relvaf2")%>%dplyr::select(Patient.ID,position,Gene,fitness),
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
  filter(size>7)%>%
  ggplot(aes(x=med.fit.plat,y=med.fit.maint,color=Gene))+
  geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.8)+
  geom_point(aes(size=size),alpha = 1)+
  geom_errorbar(aes(ymin=med.fit.maint-err.fit.maint,ymax=med.fit.maint+err.fit.maint),width=0.2)+
  geom_errorbar(aes(xmin=med.fit.plat-err.fit.plat,xmax=med.fit.plat+err.fit.plat),width=0.2)+
  scale_color_npg()+
  scale_x_continuous(name="Mean fitness Chemotherapy",lim=c(-2.5,5))+
  scale_y_continuous(name="Mean fitness maintenance",lim=c(-2.5,5))+
  scale_size(name="No. of mutations")+
  my_theme()+
  theme(legend.text = element_text(face="italic"))-> p.scatter_CvsM


png("output/figures/p.cf.scatter_CvsM.png",width=4.5, height=3,units="in",res=500,type="cairo")
p.scatter_CvsM
dev.off()

##scatterplot growthrate Platinum+PARPi vs PARPI maintenance median+iqr
full_join(df.cf_c1c7 %>% filter(variable=="relvaf2")%>%dplyr::select(Patient.ID,position,Gene,fitness),
          df.cf_c7eot%>%filter(variable=="relvaf2")%>%dplyr::select(Patient.ID,position,Gene,fitness),
          by=c("Patient.ID","Gene","position")) %>%
  filter(Patient.ID %in% (df.clin %>% filter(Number_ChemotherapyCycles > 5))$Patient.ID)%>%
  mutate(fitness_platinum = fitness.x, 
         fitness_maintenance = fitness.y)%>%
  dplyr::select(-fitness.x,-fitness.y)%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  group_by(Gene)%>%
  summarise(med.fit.plat = median(fitness_platinum),
            med.fit.maint = median(fitness_maintenance),
            min.fit.plat = quantile(fitness_platinum)[2],
            max.fit.plat = quantile(fitness_platinum)[4],
            min.fit.maint = quantile(fitness_maintenance)[2],
            max.fit.maint = quantile(fitness_maintenance)[4],
            size=n())%>%
  data.frame %>%
  filter(size>7)%>%
  ggplot(aes(x=med.fit.plat,y=med.fit.maint,color=Gene))+
  geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.8)+
  geom_point(aes(size=size))+
  geom_errorbar(aes(ymin=min.fit.maint,ymax=max.fit.maint),width=0.2)+
  geom_errorbar(aes(xmin=min.fit.plat,xmax=max.fit.plat),width=0.2)+
  scale_color_npg()+
  scale_x_continuous(name="Mean fitness Chemotherapy",lim=c(-6,6))+
  scale_y_continuous(name="Mean fitness maintenance",lim=c(-6,6))+
  theme_minimal() ->p.scatter_CvsM


png("output/figures/p.cf.scatter_CvsM.png",width=5, height=4,units="in",res=500,type="cairo")
p.scatter_CvsM
dev.off()
  
  
 