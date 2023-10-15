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

######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")
source("src/genegroup_definitions.R")




################### Mutation spectrum at day 0 ################
SC_registry %>% 
  filter(Visite == 1) %>% 
  dplyr::select(Patient.ID) %>% 
  left_join(.,df.filtered.newsamples%>%
              filter(TVAF > 0.01,tag=="true")%>%
              group_by(Patient.ID)%>%
              mutate(n.patient = n())%>%
              data.frame%>%
              dplyr::select(Patient.ID,n.patient)%>% 
              unique)%>%
  mutate(n.patient=replace_na(n.patient,0))%>%
  mutate(n.patient=ifelse(n.patient > 5, ">5",n.patient))%>%
  mutate(n.patient=factor(n.patient,levels=c("0","1","2","3","4",">5")))%>%
  dplyr::select(n.patient)%>%
  table%>%
  data.frame -> df.nom
names(df.nom) <- c("nom","Freq")

df.nom %>% 
  filter(nom!=0)%>%
  ggplot(., aes(x = nom, y = Freq, fill="1")) +
  geom_bar(stat = "identity") +
  labs(x = "No. of mutations", y = "No. of patients")+
  scale_fill_npg(name="",breaks=c(""))+
  #scale_x_continuous(limits = c(0.5, 8.5), breaks = 1:8)+
  #ggtitle("No. of mutations per patient [>1%]")+
  my_theme()->p.nom
p.nom

png("output/figures/sc/nom.png",width=4, height=3,units="in",res=500,type="cairo")
p.nom
dev.off()


##PLOTs
df.filtered.newsamples %>% 
  mutate(Gene=ifelse(Gene=="U2AF1;U2AF1L5","U2AF1",Gene))%>%
  filter(Visite==1)-> df.filtered.d0

nop <- SC_registry%>%
  filter(Visite == 1)%>%
  select(.,Patient.ID)%>%
  unique()%>%nrow #number of patients

#number of CH positive Patients
df.filtered.d0%>%
  filter(tag == "true" & TVAF >= 0.01)%>%
  select(.,Patient.ID)%>%
  unique()%>%
  nrow()

#number of mutations
df.filtered.d0 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>% 
  nrow() #number of mutations


########   Gene Mutation Prevalence Plot (plots number of gene-x-mutated patients)  #####
df.filtered.d0 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  #filter(Gene %in% ch_genes)%>%  #only CH panel, when we say: this is the prevalence plot for CH in these patients
  dplyr::select(Sample, Gene) %>% 
  data.frame %>% 
  unique %>% 
  dplyr::select(Gene) %>% 
  table %>% 
  data.frame %>% 
  filter(Freq >0) %>% 
  mutate(prev = Freq/nop) %>% 
  arrange(prev) -> prev.table
names(prev.table)<- c( "Gene","Freq","prev")

prev.table  %>%
  mutate(HRD = ifelse(is.element(Gene,hrd_genes),"HRD","non-HRD"))%>%
  ggplot(aes(x=reorder(Gene, Freq), y=prev, fill=HRD)) +
  geom_bar(stat="identity", width=0.6)+
  geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(labels = percent,limits=c(0,0.5), position = "right")+
  ylab("Gene Mutation Prevalence [%]") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  coord_flip() + 
  theme(legend.position = c(0.7, 0.3))+
  scale_fill_npg(name="",labels=c("HR-related","other")) -> p.mutprev
p.mutprev

png("output/figures/sc/mutprev.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev
dev.off()

##selected patients with multiple mutations
selected_patients <- c("SC-2","SC-3","SC-4","SC-5","SC-6","SC-7","SC-37")

df.filtered.d0 %>% 
  filter(tag == "true") %>%
  group_by(Patient.ID,Gene)%>%
  mutate(mutation = paste(Gene,row_number()))%>%
  data.frame%>%
  filter(Patient.ID %in% selected_patients) %>%
  mutate(DDR = ifelse(is.element(Gene,ddr_genes),"DDR","non-DDR"))%>%
  ggplot(aes(x=mutation, y= TVAF, fill=Gene)) +
  geom_bar(stat="identity")+
  ylab("VAF") +
  xlab("") +
  my_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())+
  facet_wrap(~Patient.ID, nrow=2)+
  scale_fill_npg() -> p.scpat
p.scpat

png("output/figures/sc/sc_patients.png",width=6, height=4,units="in",res=500,type="cairo")
p.scpat
dev.off()

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
  filter(AF < 0.01) %>%
  filter(!snp)%>%
  filter(Func=="exonic")%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.008,
         minVAF < 0.38)-> serial.mut

df.sc %>%
  mutate(Patmut = paste(Patient.ID,position,sep="_"))%>%
  filter(Patmut %in% serial.mut$Patmut)%>%
  filter(Sample_orig!="SC-3_FU")%>%  ###we don't want this in the plot
  mutate(timepoint = ifelse(Visite == 1, 0, tFU))%>%
ggplot() + 
  geom_point(aes(x=timepoint,y=TVAF,color=Gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=timepoint,y=TVAF,group=position,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=5, dir="h") +
  scale_y_log10() +
  scale_color_igv()+
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.serial


png("output/figures/sc/serial.png",width=10, height=6,units="in",res=500,type="cairo")
p.serial
dev.off()

ch_genes <- c("DNMT3A","TET2","ASXL1","PPM1D","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5")
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")

##Patient.ID
select(df.filtered,Patient.ID)%>%
  unique()->Patient.ID



