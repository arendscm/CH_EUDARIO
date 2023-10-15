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
SC_registry <- read.xlsx("data/external/SC_registry.xlsx",sheetIndex=1,header=TRUE)


# BRCA Exchange database for BRCA germline status
brcaexchange <- read.table("data/external/BRCA_Exchange_Liste_shortend.csv",sep=";",header=TRUE)

# genegroup definitions
source("src/genegroup_definitions.R")


########   Identify GERMLINE mutations-----------------------------------------------------------
##Determine BRCA status
#### BRCA Status for SC Candidates ####
df->df.orig
df.orig %>% filter(is.na(Patient.ID)) %>% dplyr::select(-Patient.ID)%>%full_join(.,SC_registry %>% dplyr::select(Sample_orig,Patient.ID),by=c("Sample_orig")) -> df

df %>% 
  filter(Gene == "BRCA1"|Gene =="BRCA2") %>%
  filter(TVAF > 0.25) %>%
  filter(ExonicFunc!= "synonymous SNV") %>%
  filter(AF < 0.05)%>%
  left_join(., brcaexchange, by="Genomic_Coordinate_hg38") %>%
  filter(is.element(BRCA.Exchange_Pathogenicity_expert,c("Pathogenic","Not Yet Reviewed"))|
           (is.na(BRCA.Exchange_Pathogenicity_expert)&
              (is.element(ExonicFunc,c("frameshift substitution","stopgain","startloss"))|
                 is.element(Func,c("splicing","exonic;splicing")))))%>%## all variants that are classified as pathogenic by expert panel or not yet reviewed, or that have no match in BRCA exchange but are truncating
  filter(!str_detect(BRCA.Exchange_Clinical_Significance_ClinVar,"Benign")|is.na(BRCA.Exchange_Clinical_Significance_ClinVar))-> df.brca_germline_SC


SC_registry %>% 
  filter(Visite==1) %>%
  mutate(brca1_germline = ifelse(is.element(Patient.ID,(df.brca_germline_SC%>%filter(Gene=="BRCA1"))$Patient.ID),1,0))%>%
  mutate(brca2_germline = ifelse(is.element(Patient.ID,(df.brca_germline_SC%>%filter(Gene=="BRCA2"))$Patient.ID),1,0))%>%
  dplyr::select(Patient.ID,brca1_germline,brca2_germline) %>% 
  unique -> id.brca_germline_SC

##identify other HRD Gene germline mutations

df %>% 
  filter(is.element(Gene,hrd_genes)) %>%
  filter(TVAF > 0.25) %>%
  filter(ExonicFunc!= "synonymous SNV") %>%
  filter(Func == "exonic")%>%
  filter(AF < 0.01) %>%
  filter(!is.element(Gene,c("BRCA1","BRCA2")))%>%
  filter(is.element(ExonicFunc,c("frameshift substitution","stopgain","startloss"))) -> df.hrd_germline_SC # list has to be discussed with an expert

SC_registry %>% 
  filter(Visite==1) %>%
  mutate(hrd_germline = ifelse(is.element(Patient.ID,df.hrd_germline_SC$Patient.ID),1,0))%>%
  dplyr::select(Patient.ID,hrd_germline) %>% 
  unique -> id.hrd_germline_SC


tempdata <-ls()
rm(list=tempdata[!is.element(tempdata,c("df.brca_germline_SC","id.brca_germline_SC","df.hrd_germline_SC","id.hrd_germline_SC"))])
rm(tempdata)

save.image("data/interim/brca_sc.RData")

#filename="output/SC-BRCAandHRD-Status.xlsx"
#write.xlsx(df.brca_germline_SC,filename,sheetName="BRCA",append=TRUE)
#write.xlsx(df.hrd_germline_SC,filename, sheetName="HRD",append=TRUE)



