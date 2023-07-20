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
library(ggsci)

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')
load('data/interim/brca.RData')

######## Global input
source("src/material_table.R") #Patient ids
source("src/genegroup_definitions.R") #gene group definitions
source("src/global_functions_themes.R") #functions and themes

########   SERIAL SAMPLES

#### identity check via SNPs ####
df %>% 
  filter(Visite=="EOT") %>% 
  filter(Material=="wb")%>%
  mutate(ID=paste(Patient.ID,position))-> df.eot1

ids %>% filter(Visite=="EOT",Material == "wb") -> eotsamples

df %>% 
  filter(Visite=="C1D1", Material == "wb") %>% 
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
  filter(AF<0.01)%>%
  filter(snp==FALSE)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         minp.binom = min(p.binom)) %>%
  data.frame()%>%
  filter(maxVAF > 0.008) %>%
  filter(minVAF < 0.4)%>%
  filter(minp.binom ==-Inf) %>%
  filter(Patient.ID != "4221006")%>% ##very short FU time
  .$Patmut ->Patmut.serial2  #some "partners" get kicked out here, rescue them back

df.eot %>% filter(tag == "true") %>% .$Patmut -> Patmut.serial1 ##include also mutations called in just one timepoint, but manually checked in igv

Patmut.serial <- c(Patmut.serial1,Patmut.serial2)%>% unique

##some mutations are only called in one timepoint and were beneath the threshold of 0.001 VAF in the other. 
##Here we include them and set their VAF to the detection threshold

df.eot %>%
  filter(Patmut %in% Patmut.serial) %>%
  group_by(Patmut) %>%
  mutate(n.serial = n())%>%
  data.frame%>%
  filter(n.serial == 1) %>%
  mutate(Visite = ifelse(Visite=="C1D1","EOT","C1D1"),
         TVAF = 0.001,
         timepoint = ifelse(Visite=="C1D1",0,tEOT)) %>%
  dplyr::select(-n.serial) -> df.missing

df.serial <- full_join(df.eot,df.missing)
  
df.serial%>%
  filter(Patmut %in% Patmut.serial)%>%
  ggplot() + 
  geom_point(aes(x=timepoint,y=TVAF,color=Gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=TVAF,group=position,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=9, dir="h") +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_y_log10()+
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  scale_color_igv()+
  theme_minimal()-> p.serial
p.serial

png("output/figures/p.serial_wb.png",width=12, height=8,units="in",res=500,type="cairo")
p.serial
dev.off()

##same plot with CH genes only

df.serial%>%
  filter(Gene %in% typical_ch_genes)%>%
  filter(Patmut %in% Patmut.serial)%>%
  ggplot() + 
  geom_point(aes(x=timepoint,y=TVAF,color=Gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=TVAF,group=position,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=8, dir="h") +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_y_log10()+
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  scale_color_npg()+
  theme_minimal()-> p.serial_ch
p.serial_ch

png("output/figures/p.serial_wb_ch.png",width=12, height=8,units="in",res=500,type="cairo")
p.serial_ch
dev.off()


#### SERIAL SAMPLES growth preprocessing ####
df.serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "C1D1") %>%
  mutate(vaf_d1=TVAF,
         depth_d1=readDepth,
         vc_d1=TR2) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,depth_d1,vc_d1,vaf_d1)-> df.eotd1

df.serial %>% 
  filter(Patmut %in% Patmut.serial)%>%
  filter(Visite == "EOT") %>%
  mutate(vaf_eot=TVAF,
         depth_eot=readDepth,
         vc_eot=TR2) %>%
  dplyr::select(Patient.ID,Gene,AAChange,ExonicFunc,position,depth_eot,vc_eot,vaf_eot,timepoint)-> df.eoteot

#function to calculate fitness error
h <- function(x){
  return(1/x +(2/(1-2*x)))
}

df.eot_rel <- full_join(df.eotd1,df.eoteot) %>% mutate(relvaf1 = vaf_d1/vaf_d1,
                                                       relvaf2 = vaf_eot/vaf_d1) %>%
  mutate(dt=timepoint/365,
         sig.x = sqrt((vaf_d1*(1-vaf_d1))/depth_d1), 
         sig.y = sqrt((vaf_eot*(1-vaf_eot))/depth_eot),
         fitness = log(vaf_eot/vaf_d1*(1-2*vaf_d1)/(1-2*vaf_eot))/dt,
         dfitness = 1/dt*sqrt((h(vaf_d1)*sig.x)^2+(h(vaf_eot)*sig.x)^2),
         dyn = ifelse((fitness-1.96*dfitness)>0,1,ifelse((fitness+1.96*dfitness)<0,-1,0)))%>%
  melt.data.frame(measure.vars = c("relvaf1","relvaf2"))


####   plot rel vaf2 as points according to gene, coloured in ExonicFunc (frameshift,...) ####
df.eot_rel %>% 
  mutate(timepoint = ifelse(variable=="relvaf1",0,timepoint))%>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  mutate(fitness_binom = ifelse(fitness>0.25,"increasing",
                                ifelse(fitness < -0.25,"decreasing","stable")))%>%
  filter(is.element(Gene,c("CHEK2","PPM1D","DNMT3A","TP53","TET2", "ATM")))%>%
  ggplot(aes(x = timepoint, y = value, color = fitness_binom, group = Patmut)) +
  geom_point(size = 1.5, na.rm = FALSE) + 
  geom_line()+
  scale_y_log10() +
  scale_color_manual(values=c("#4DBBD5FF","#E64B35FF","#00A087FF"))+
  labs(x = "Time in days", y = "VAF change", color = "Fitness") +
  my_theme() +
  theme(strip.text = element_text(face = "italic"))+
  facet_wrap(~Gene)-> p.growth
p.growth

png("output/figures/relgrowth_wb_gene.png",width=8, height=5,units="in",res=500,type="cairo")
p.growth
dev.off()

####   Boxplot Fitness index according to DDR/non DDR ####
my_comp=list(c("DNMT3A","TP53"),c("PPM1D","TET2"),c("DNMT3A","PPM1D"))

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
            ylim = c(-4,16),
            size=0.8,
            alpha=1,
            repel=TRUE,
            #yscale = "log10",
            scales = "free",
            add = c("jitter")
           )+
  stat_compare_means(comparisons=my_comp,label="p.signif",vjust=0.01,label.y = c(10,12,14))+
  scale_color_npg()+
  theme_minimal() + 
  theme(axis.title.x = element_blank()) +
  theme(#legend.position = "none",
        axis.text.x = element_text(face="italic"),
        axis.title.y = element_text(face ="plain"),
        plot.title = element_text(hjust=0,face ="plain")) ->p.fitness_boxplot

png("output/figures/fitness_boxplot.png",width=5, height=5,units="in",res=500,type="cairo")
p.fitness_boxplot
dev.off()

##alternative: violin plot
my_comp=list(c("DNMT3A","TP53"),c("DNMT3A","PPM1D"),c("TET2","TP53"),c("TET2","PPM1D"))

df.eot_rel %>% 
  filter(variable == "relvaf2") %>% 
  group_by(Gene)%>%
  mutate(medfit = median(fitness))%>%
  data.frame%>%
  filter(is.element(Gene,c("CHEK2","PPM1D","DNMT3A","TP53","TET2")))%>%
  mutate(DDR = ifelse(is.element(Gene,c("TP53","PPM1D","CHEK2")),"DDR","non-DDR"))%>%
  ggplot(., aes(x=reorder(Gene,medfit), y=fitness, fill=DDR)) + 
  geom_violin(trim=TRUE,bw=1,aes(color=DDR))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="Gene", y = "Fitness")+
  my_theme() + 
  theme(axis.title.x = element_blank()) +
  ylim(c(-5,16))+
  scale_fill_npg(name="",labels=c("DDR gene","DTA gene"))+
  scale_color_npg(name="",labels=c("DDR gene","DTA gene"))+
  geom_hline(yintercept=0, linetype="dashed",alpha=0.7)+
  theme(#legend.position = "none",
    axis.text.x = element_text(face="italic"),
    axis.title.y = element_text(face ="plain"),
    plot.title = element_text(hjust=0,face ="plain"))+
  stat_compare_means(comparisons=my_comp,
                     label="p.signif",
                     vjust=0.001,
                     label.y = c(10,11,12,13),
                     tip.length=c(0.005,0.005,0.005,0.005))->p.fitness_violin

png("output/figures/fitness_violin.png",width=6, height=6,units="in",res=500,type="cairo")
p.fitness_violin
dev.off()


##emerging and disappearing clones crossing the 1% VAF threshold
df.eot_rel %>% filter(vaf_d1 < 0.01&vaf_eot >= 0.01) %>% mutate(class = "emerging")-> df.emerging
df.eot_rel %>% filter(vaf_d1 >= 0.01&vaf_eot < 0.01) %>% mutate(class = "disappearing")-> df.disappearing

full_join(df.emerging,df.disappearing) %>%
  mutate(gene = ifelse(is.element(Gene,hrd_genes),"HRD",Gene))%>%
  ggplot(aes(x=class,fill=gene))+
  geom_bar(stat="count",width=0.5)+
  labs(title="",x="", y = "No. of mutations")+
  my_theme() + 
  theme(axis.title.x = element_blank(),
        legend.text = element_text(face="italic")) +
  scale_fill_npg(name="Gene") -> p.dynamics
  
png("output/figures/barchart_dynamics.png",width=4, height=4,units="in",res=500,type="cairo")
p.dynamics
dev.off()

#### Fitness of PPM1D by brca status (to see whether dynamics unter PARP Inhb. differ depending on BRCA status?) #### 
df.eot_rel %>% 
  filter(variable == "relvaf2") %>% 
  group_by(Gene)%>%
  mutate(medfit = median(fitness))%>%
  data.frame%>%
  filter(is.element(Gene,c("PPM1D")))%>%
  left_join(.,id.brca_germline,by = "Patient.ID")%>%
  mutate(brca_germline=brca1_germline+brca2_germline)%>%
  ggplot(., aes(x=factor(brca_germline), y=fitness, group=as.factor(brca_germline))) + 
  geom_violin(trim=TRUE,bw=1)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="Gene", y = "Fitness")+
  my_theme() + 
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed",alpha=0.7)+
  theme(#legend.position = "none",
    axis.title.y = element_text(face ="plain"),
    plot.title = element_text(hjust=0,face ="plain"))+
  stat_compare_means()


