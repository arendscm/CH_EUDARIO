# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: mutational spectrum in patients with and without prior PARPi
#
# Input: preprocessed clinical data and mutational data
#
# Output: data.frame with clinical and genomic data
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
library(dplyr)
library(tableone)
library(ggplot2)
library(ggsci)

########   Load IDs    ########
source("src/material_table.R")
source("src/global_functions_themes.R")

########## Load clinical data  ##########
load("data/interim/clin.RData")
load("data/interim/seqdata_filtered.RData")

######## Plot with mut spectrum in pts with and without prior PARPi

df.clin %>% mutate(BRCA_germline = brca1_germline+brca2_germline) -> df.clin
df.clin$BRCA_germline %>% table %>% data.frame -> n.brca
names(n.brca) <- c("BRCA_germline","n.brca")

genes=c("DNMT3A","TET2","ASXL1","PPM1D","TP53","CHEK2","ATM")

df.filtered.c1d1 %>% left_join(.,df.clin %>% 
                                 group_by(BRCA_germline) %>% 
                                 mutate(n.brca = n()) %>% 
                                 data.frame %>% 
                                 dplyr::select(Patient.ID, BRCA_germline,n.brca))%>%
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  filter(Gene %in% genes)%>%
  dplyr::select(Sample, Gene, BRCA_germline,n.brca) %>% 
  data.frame %>% 
  unique %>% 
  dplyr::select(Gene,BRCA_germline) %>% 
  table %>% 
  data.frame %>% 
  left_join(.,n.brca)%>%
  mutate(prev = Freq/n.brca) %>% 
  arrange(prev) %>%
  ggplot(aes(x=reorder(Gene, Freq), y=prev, fill=BRCA_germline)) +
  geom_bar(stat="identity", width=0.6, position=position_dodge())+
  xlab("")+
  scale_y_continuous(labels = percent,limits=c(0,0.4), position = "right")+
  ylab("Gene Mutation Prevalence [%]") +
  my_theme()+
  scale_fill_npg()+
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  coord_flip() -> p.mutprev_brca
p.mutprev_brca

png("output/figures/mutprev_brca.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev_brca
dev.off()



