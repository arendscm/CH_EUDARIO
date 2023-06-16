# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: Defines BRCA Germline Mutations (BRCA/HRD and their IDs)
#
# Input: seqdata
#
# Output: df.brca_germline and BRCA germline IDs/ df.hrd_germline
#
# ==============================================================================
########   Dependencies   #####
library(dplyr)
library(base)
library(reshape)
library(tidyr)
library(stringr)
library(xlsx)

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')

######## Get sc registry
##

# BRCA Exchange database for BRCA germline status
brcaexchange <- read.table("data/external/BRCA_Exchange_Liste_shortend.csv",sep=";",header=TRUE)

########   Identify GERMLINE mutations-----------------------------------------------------------
##Determine BRCA status
#### BRCA Status for SC Candidates ####
df->df.orig
df.orig %>% filter(is.na(Patient.ID)) -> df

df %>% 
  filter(Gene == "BRCA1"|Gene =="BRCA2") %>%
  filter(TVAF > 0.25) %>%
  filter(ExonicFunc!= "synonymous SNV") %>%
  filter(AF < 0.05)%>%
  left_join(., brcaexchange, by="Genomic_Coordinate_hg38") %>%
  filter(is.element(BRCA.Exchange_Pathogenicity_expert,c("Pathogenic","Not Yet Reviewed"))|
           (is.na(BRCA.Exchange_Pathogenicity_expert)&
              (is.element(ExonicFunc,c("frameshift substitution","stopgain"))|
                 is.element(Func,c("splicing","exonic;splicing")))))%>%## all variants that are classified as pathogenic by expert panel or not yet reviewed, or that have no match in BRCA exchange but are truncating
  filter(!str_detect(BRCA.Exchange_Clinical_Significance_ClinVar,"Benign")|is.na(BRCA.Exchange_Clinical_Significance_ClinVar))-> df.brca_germline_SC

##identify other HRD Gene germline mutations
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")

df %>% 
  filter(is.element(Gene,hrd_genes)) %>%
  filter(TVAF > 0.25) %>%
  filter(ExonicFunc!= "synonymous SNV") %>%
  filter(Func == "exonic")%>%
  filter(AF < 0.01) %>%
  filter(is.element(ExonicFunc,c("frameshift substitution","stopgain"))) -> df.hrd_germline_SC # list has to be discussed with an expert

filename="output/SC-BRCAandHRD-Status.xlsx"
write.xlsx(df.brca_germline_SC,filename,sheetName="BRCA",append=TRUE)
write.xlsx(df.hrd_germline_SC,filename, sheetName="HRD",append=TRUE)



