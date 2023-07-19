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

#load brca exchange database
brcaexchange <- read.table("data/external/BRCA_Exchange_Liste_shortend.csv",sep=";",header=TRUE)


#####BRCA germline and somatic mutations c1d1
df.c1d1 %>%
  filter(Gene == "BRCA1"|Gene == "BRCA2") %>%
  filter(ExonicFunc != "synonymous SNV")%>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 19) %>%
  filter(TVAF >=0.01) %>%
  filter(p.binom <= -Inf)%>%
  filter(mutFreq < 0.1*n.lane)%>%
  filter(AF < 0.05) %>%
  filter(!snp) %>% 
  left_join(., brcaexchange, by="Genomic_Coordinate_hg38") %>%
  filter(is.element(BRCA.Exchange_Pathogenicity_expert,c("Pathogenic","Not Yet Reviewed"))|
           (is.na(BRCA.Exchange_Pathogenicity_expert)&
              (is.element(ExonicFunc,c("frameshift substitution","stopgain","startloss"))|
                 is.element(Func,c("splicing","exonic;splicing")))))%>%## all variants that are classified as pathogenic by expert panel or not yet reviewed, or that have no match in BRCA exchange but are truncating
  filter(!str_detect(BRCA.Exchange_Clinical_Significance_ClinVar,"Benign")|is.na(BRCA.Exchange_Clinical_Significance_ClinVar))%>%
  filter(TVAF < 0.2)->df.brca_somatic

ids %>% 
  filter(Material=="cf") %>%
  filter(Visite == "C1D1")%>%
  mutate(brca1_somatic = ifelse(is.element(Patient.ID,(df.brca_somatic%>%filter(Gene=="BRCA1"))$Patient.ID),1,0))%>%
  mutate(brca2_somatic = ifelse(is.element(Patient.ID,(df.brca_somatic%>%filter(Gene=="BRCA2"))$Patient.ID),1,0))%>%
  dplyr::select(Patient.ID,brca1_somatic,brca2_somatic) %>% 
  unique -> id.brca_somatic

tempdata <-ls()
rm(list=tempdata[tempdata != "id.brca_somatic"&tempdata != "df.brca_somatic"])
rm(tempdata)

save.image("data/interim/brca_somatic.RData")