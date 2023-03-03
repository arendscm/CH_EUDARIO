# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: from variantcalls to excel list of filtered results
#
# Input: seqdata
#
# Output: Excel list of filtered results
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

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########  Load preprocessed sequencing data
load('data/interim/seqdata.RData')

##load hotspot data
load('data/interim/hotspots.RData')

######## Get Patient ids
source("src/ids.R")

########   FILTERING cf calls------------------------------------------------------------
df -> df.backup
df.backup %>% filter(Material == "cf",
              Visite == "C1D1") -> df.c1d1

df.backup %>% filter(Material == "cf",
                 Visite == "EOT") -> df.eot

##HRD genes
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3", "BRCA1", "BRCA2")

##load BRCA germline mutations
source("src/brca_germline.R")

#####BRCA germline and somatic mutations c1d1
df.c1d1 %>%
  filter(Gene == "BRCA1"|Gene == "BRCA2") %>%
  #filter(ExonicFunc != "synonymous SNV")%>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 19) %>%
  filter(TVAF >0.01) %>%
  filter(p.binom< -10)%>%
  filter(mutFreq < 0.1*n.lane)%>%
  filter(!snp)

##brca germline and somatic mutations at EOT
df.eot %>%
  filter(Gene == "BRCA1"|Gene == "BRCA2") %>%
  #filter(ExonicFunc != "synonymous SNV")%>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 19) %>%
  filter(TVAF >0.005) %>%
  filter(p.binom< -10)%>%
  filter(mutFreq < 0.1*n.lane)%>%
  filter(!snp)%>%select(Patient.ID,Chr,Start,End,Ref,Alt,Gene,Func,ExonicFunc,TVAF,n.material,n.visite,c1d1_cf,eot_cf)
