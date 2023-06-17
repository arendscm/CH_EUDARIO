# ______________________________________________________________________________
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: Analysis of longitudinal data WB samples
#
# Input: seqdata
#
# Output: serial sample analysis of cf DNA and wb DNA mutations
# 
# ______________________________________________________________________________
########   Dependencies   #####
library(base)
library(dplyr)
library(stringr)
library(reshape)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(reshape)
library(ggpubr)
library(maftools)

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')

######## Global input
source("src/material_table.R") #Patient ids
source("src/genegroup_definitions.R") #gene group definitions
source("src/global_functions_themes.R") #functions and themes

########   SERIAL SAMPLES

#### identity check via SNPs ####
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

#### serial samples preparation ####
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

df %>% 
  filter(!is.na(Patient.ID))%>%
  filter(is.na(replicate))%>%
  filter(c1d1_wb==1&eot_wb==1) %>% 
  filter(Material=="wb")%>%
  mutate(timepoint = ifelse(Visite == "C1D1",0,ifelse(Visite == "EOT",tEOT,NA)))-> df.eot

#### serial samples dynamics plot ####
df.eot %>% 
  filter(n.material>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(p.binom < -10) %>% 
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.008) %>%
  filter(minVAF < 0.35)%>%
  .$Patmut ->Patmut.serial2  #some "partners" get kicked out here, rescue them back

df.eot%>%
  filter(Patmut %in% Patmut.serial2)%>%
  ggplot() + 
  geom_point(aes(x=timepoint,y=TVAF,color=Gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=TVAF,group=position,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=10, dir="h") +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_y_log10()+
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  scale_color_igv()+
  theme_minimal()-> p.serial

png("output/figures/p.serial.png",width=8, height=2,units="in",res=500,type="cairo")
p.serial
dev.off()

##Examples
df.eot %>% 
  filter(Patmut %in% Patmut.serial2)%>%
  filter(is.element(Patient.ID,c("4202011","4207004","4202001","4221010")))%>%
  ggplot() + 
  geom_point(aes(x=timepoint,y=TVAF,color=Gene,group=Patient.ID),size=1 ,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=TVAF,group=position,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=6, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.45)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.serial

png("output/figures/p.serial_example.png",width=6, height=3,units="in",res=500,type="cairo")
p.serial
dev.off()


#### Serial samples by brca status (question: do dynamics unter PARP Inhb. differ depending on BRCA status?) ####
df.eot %>% 
  filter(Patmut %in% Patmut.serial2)%>%
  left_join(.,id.brca_germline,by = "Patient.ID")%>%
  mutate(brca_germline=brca1_germline+brca2_germline)%>%
  filter(Gene=="PPM1D"|Gene=="TP53")%>%
  ggplot() + 
  geom_point(aes(x=timepoint,y=TVAF,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=TVAF,group=position,color=Gene),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ brca_germline, ncol=2, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.serial_brca

#### SERIAL SAMPLES growth preprocessing ####
df.eot %>% 
  filter(Patmut %in% Patmut.serial2)%>%
  filter(Visite == "C1D1") %>%
  mutate(vaf_d1=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_d1)-> df.eotd1

df.eot %>% 
  filter(Patmut %in% Patmut.serial2)%>%
  filter(Visite == "EOT") %>%
  mutate(vaf_eot=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,vaf_eot,timepoint)-> df.eoteot

df.eot_rel <- full_join(df.eotd1,df.eoteot) %>% mutate(relvaf1 = vaf_d1/vaf_d1,
                                                       relvaf2 = vaf_eot/vaf_d1) %>%
  mutate(fitness = log(vaf_eot/vaf_d1*(1-2*vaf_d1)/(1-2*vaf_eot))/(timepoint/365))%>%
  melt.data.frame(measure.vars = c("relvaf1","relvaf2"))


####   growth when first vaf set to value 1.0  ####
df.eot_rel %>%  
  ggplot() + 
  geom_point(aes(x=timepoint,y=value,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=value,group=position,color=Gene),size=1*1,na.rm=FALSE) + 
  scale_y_log10() +
  labs(x="Timepoint",y="log(VAF change)",colour="Mutated Gene") +
  theme_minimal()-> p.vafchange

####   plot rel vaf2 as points according to gene ####
df.eot_rel %>% 
  filter(variable == "relvaf2") %>% 
  ggplot(aes(x = reorder(Gene, value, FUN = median), y = value, color = Gene, group = Patient.ID)) +
  geom_point(size = 1.5, na.rm = FALSE) + 
  scale_y_log10() +
  labs(x = "Gene", y = "log(VAF change)", colour = "Mutated Gene") +
  theme_graphicalabstract() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))-> p.growth

png("output/figures/p.relvaf2.png",width=6, height=4,units="in",res=500,type="cairo")
p.serial
dev.off()

####   plot rel vaf2 as points according to gene, coloured in ExonicFunc (frameshift,...) ####
df.eot_rel %>% 
  filter(is.element(Gene,c("CHEK2","PPM1D","DNMT3A","TP53","TET2", "ATM")))%>%
  filter(variable == "relvaf2") %>% 
  ggplot(aes(x = reorder(Gene, value, FUN = median), y = value, color = ExonicFunc, group = Patient.ID)) +
  geom_point(size = 1.5, na.rm = FALSE) + 
  scale_y_log10() +
  labs(x = "Gene", y = "log(VAF change)", colour = "Mutated Gene") +
  theme_graphicalabstract() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))-> p.growth

png("output/figures/p.relvaf_ExonicFunc.png",width=6, height=4,units="in",res=500,type="cairo")
p.serial
dev.off()

####   Boxplot Fitness index according to DDR/non DDR ####
my_comp=list(c("DNMT3A","PPM1D"),c("DNMT3A","TP53"),c("DNMT3A","CHEK2"),c("PPM1D","TET2"),c("TET2","TP53"))

df.eot_rel %>% 
  filter(variable == "relvaf2") %>% 
  #mutate(growthrate = log(value)/timepoint)%>%
  filter(is.element(Gene,c("CHEK2","PPM1D","DNMT3A","TP53","TET2")))%>%
  mutate(DDR = ifelse(is.element(Gene,c("TP53","PPM1D","CHEK2")),"DDR","non-DDR"))%>%
  ggboxplot(., 
            x = "Gene",
            y = "fitness",
            order = c("TP53","PPM1D","CHEK2","TET2","DNMT3A"),
            combine = TRUE,
            color = "DDR", 
           # palette = ,
            xlab = "Gene",
            ylab = "Fitness",
            title = "",
            width = 0.3,
            #ylim = c(-0.01,0.02),
            size=0.8,
            alpha=1,
            repel=TRUE,
            #yscale = "log10",
            scales = "free",
            add = c("jitter")
           )+
  stat_compare_means(comparisons=my_comp)+
  theme_minimal() + 
  theme(axis.title.x = element_blank()) +
  theme(#legend.position = "none",
        axis.text.x = element_text(face="italic"),
        axis.title.y = element_text(face ="plain"),
        plot.title = element_text(hjust=0,face ="plain")) ->p.growth

png("output/figures/growth.png",width=6, height=4,units="in",res=500,type="cairo")
p.growth
dev.off()

##alternative: violin plot
df.eot_rel %>% 
  filter(variable == "relvaf2") %>% 
  group_by(Gene)%>%
  mutate(medfit = median(fitness))%>%
  data.frame%>%
  filter(is.element(Gene,c("CHEK2","PPM1D","DNMT3A","TP53","TET2")))%>%
  mutate(DDR = ifelse(is.element(Gene,c("TP53","PPM1D","CHEK2")),"DDR","non-DDR"))%>%
  ggplot(., aes(x=reorder(Gene,medfit), y=fitness, fill=DDR)) + 
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="Gene", y = "Fitness")+
  my_theme() + 
  theme(axis.title.x = element_blank()) +
  ylim(c(-3,8))+
  scale_fill_npg()+
  theme(#legend.position = "none",
    axis.text.x = element_text(face="italic"),
    axis.title.y = element_text(face ="plain"),
    plot.title = element_text(hjust=0,face ="plain")) ->p.growth

png("output/figures/growth.png",width=6, height=4,units="in",res=500,type="cairo")
p.growth
dev.off()



