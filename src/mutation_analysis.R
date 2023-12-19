# ______________________________________________________________________________
# CH in EUDARIO
#
# Author: Max & Klara
#
# Description: Mutation analysis and plots from filtered mutational data
#
# Input: df.filtered from data/interim/seqdata_filtered.RData
#
# Output: Excel list of filtered results, plots, ...
# 
# ______________________________________________________________________________
########   Dependencies and load data  ####
library(base)
library(dplyr)
library(stringr)
library(reshape)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(reshape)
library(ggpubr)
library(g3viz)
library(maftools)
library(RColorBrewer)
library(ggsci)

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata_filtered.RData')

######## Get Patient ids
source("src/material_table.R")

######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")
source("src/genegroup_definitions.R")



######## Available samples Overview  ##############

df.material %>% 
  dplyr::select(-sum_wb,-sum_cf)%>%
  mutate(eot_ue_cf=ue_cf+eot_cf)%>%
  dplyr::select(-eot_cf,-ue_cf)%>%
  melt.data.frame(id="Patient.ID") %>%
  ggplot(aes(y=variable,x=Patient.ID,fill=as.factor(value)))+
  geom_tile()+
  scale_fill_npg(name="Sample sequenced")+
  xlab("Visit")+
  ylab("Patient ID")+
  my_theme()+
  theme(axis.text.x=element_text(angle=90))->p.material

png("output/figures/material.png",width=15, height=7,units="in",res=500,type="cairo")
p.material
dev.off()

#get numbers
df.material %>% mutate(eot_ue_cf = eot_cf+ue_cf) %>% dplyr::select(-Patient.ID) %>% sign() %>% colSums()

########   Mutational Analysis preparation #####

##PLOTs
df.filtered.c1d1 %>% mutate(Gene=ifelse(Gene=="U2AF1;U2AF1L5","U2AF1",Gene)) -> df.filtered.c1d1

nop <- ids%>%
  filter(Visite == "C1D1" & Material == "wb")%>%
  select(.,Patient.ID)%>%
  unique()%>%nrow #number of patients

#number of mutated positive Patients
df.filtered.c1d1%>%
  filter(tag == "true" & TVAF >= 0.01)%>%
  select(.,Patient.ID)%>%
  unique()%>%
  nrow()

#number of mutations
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>% 
  nrow() #number of mutations

#number of typical CH mutations
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  filter(Gene %in% typical_ch_genes)%>%
  nrow()

#number of HRD mutations
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  filter(Gene %in% hrd_genes)%>%
  nrow()

#number of HRD + non-CH mutations
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  filter(Gene %in% setdiff(hrd_genes,ch_genes))%>%
  nrow()

########   Gene Mutation Prevalence Plot (plots number of gene-x-mutated patients)  #####
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  mutate(gene_group = ifelse(Gene %in% typical_ch_genes_without_HRD,"CH",
                             ifelse(Gene %in% typical_ch_genes,"CH / HR-related",
                                    ifelse(Gene %in% hrd_genes,"HR-related","other myeloid")))) %>%
  dplyr::select(Sample, Gene, gene_group) %>% 
  data.frame %>% 
  unique %>% 
  dplyr::select(Gene, gene_group) %>% 
  table %>% 
  data.frame %>% 
  filter(Freq >0) %>% 
  mutate(prev = Freq/nop) %>% 
  arrange(prev) -> prev.table
names(prev.table)<- c( "Gene","gene_group","Freq","prev")

prev.table  %>%
  ggplot(aes(x=reorder(Gene, -Freq), y=prev, fill=gene_group)) +
  geom_bar(stat="identity", width=0.7)+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(labels = percent,limits=c(0,0.3))+
  ylab("Gene mutation prevalence") +
  my_theme() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.35,face="italic"),
        axis.ticks.x = element_blank()) +
  theme(legend.position = c(0.8, 0.8))+
  scale_fill_manual(values = c("#E64B35FF","#3C5488FF","#4DBBD5FF","#00A087FF"), name="Gene group") -> p.mutprev
p.mutprev

png("output/figures/mutprev2.png",width=8, height=4,units="in",res=500,type="cairo")
p.mutprev
dev.off()


########   Number of Mutations per Gene with type of mutation  ####

#########    No of mutations per patient
ids %>% 
  filter(visit_material=="C1D1_wb") %>% 
  dplyr::select(Patient.ID) %>% 
  left_join(.,df.filtered.c1d1%>%
              filter(TVAF >= 0.01,tag=="true",Gene %in% typical_ch_genes)%>%
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

png("output/figures/nom.png",width=3, height=3,units="in",res=500,type="cairo")
p.nom
dev.off()


########   mutation frequency according to type of mutation  ####
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%  
  filter(TVAF >= 0.01) %>%
  dplyr::select(Gene,ExonicFunc) %>% 
  mutate(ExonicFunc = replace(ExonicFunc,ExonicFunc == ".","splice mutation")) %>%
  data.frame %>% 
  group_by(Gene) %>% 
  summarise(Gene.freq = n())%>%
  as.data.frame %>% 
  full_join(df.filtered.c1d1 %>% filter(tag == "true"),.,by = "Gene") -> df.genefreq

df.genefreq %>%
  dplyr::select(Gene,ExonicFunc) %>% 
  mutate(ExonicFunc = replace(ExonicFunc,ExonicFunc == ".","splice mutation")) %>%
  data.frame %>% 
  table %>% 
  data.frame %>%
  ggplot(aes(x=reorder(Gene,Freq), y=Freq,fill=ExonicFunc)) +
  geom_bar(stat="identity", width=0.6)+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(position = "right")+
  ylab("Mutation Frequency") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  scale_fill_npg() +
  labs(fill="")+
  coord_flip()+
  theme(legend.position = c(0.7, 0.3))-> p.mutfreq
p.mutfreq


png("output/figures/mutfreq.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutfreq
dev.off()

##create xlsx for supplements
variables <- c("Patient.ID","Chr","Start","End","Ref","Alt","Gene","Func","ExonicFunc","AAChange","cosmic92_coding","readDepth","TR1","TR2","TVAF")
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  mutate(gene_group = ifelse(Gene %in% typical_ch_genes_without_HRD,"CH",
                             ifelse(Gene %in% typical_ch_genes,"CH / HR-related",
                                    ifelse(Gene %in% hrd_genes,"HR-related","other myeloid")))) %>%
  dplyr::select(variables) -> df.mut

write.xlsx(df.mut, "output/tables/mutations_d1.xlsx",sheetName="SomaticMutations_d1")
