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

########  Load preprocessed sequencing data
load('data/interim/seqdata.RData')

######## Get Patient ids
source("src/ids.R")
source("src/genegroup_definitions.R")

# BRCA Exchange database for BRCA germline status
brcaexchange <- read.table("data/external/BRCA_Exchange_Liste_shortend.csv",sep=";",header=TRUE)

########   Identify GERMLINE mutations-----------------------------------------------------------
##Determine BRCA status
df %>% 
  filter(!is.na(Patient.ID))%>%
  filter(is.na(replicate))%>%
  filter(firstTimepoint_wb == 1)%>%
  filter(Material == "wb")%>%
  filter(Gene == "BRCA1"|Gene =="BRCA2") %>%
  filter(TVAF > 0.1) %>%
  filter(ExonicFunc!= "synonymous SNV") %>%
  filter(tag!="false")%>%
  filter(AF < 0.05)%>%
  left_join(., brcaexchange, by="Genomic_Coordinate_hg38") %>%
  filter(is.element(BRCA.Exchange_Pathogenicity_expert,c("Pathogenic","Not Yet Reviewed"))|
           (is.na(BRCA.Exchange_Pathogenicity_expert)&
              (is.element(ExonicFunc,c("frameshift substitution","stopgain","startloss"))|
                 is.element(Func,c("splicing","exonic;splicing")))))%>%## all variants that are classified as pathogenic by expert panel or not yet reviewed, or that have no match in BRCA exchange but are truncating
  filter(!str_detect(BRCA.Exchange_Clinical_Significance_ClinVar,"Benign")|is.na(BRCA.Exchange_Clinical_Significance_ClinVar))-> df.brca_germline


##create list of all patients incl BRCA status
ids %>% 
  filter(Material=="wb") %>%
  filter(firstTimepoint_wb == 1)%>%
  mutate(brca1_germline = ifelse(is.element(Patient.ID,(df.brca_germline%>%filter(Gene=="BRCA1"))$Patient.ID),1,0))%>%
  mutate(brca2_germline = ifelse(is.element(Patient.ID,(df.brca_germline%>%filter(Gene=="BRCA2"))$Patient.ID),1,0))%>%
  dplyr::select(Patient.ID,brca1_germline,brca2_germline) %>% 
  unique -> id.brca_germline ##this is a list of patient ids with brca status

##identify other HRD Gene germline mutations

df %>% 
  filter(!is.na(Patient.ID))%>%
  filter(Visite == "C1D1")%>%
  filter(Material == "wb")%>%
  filter(is.element(Gene,hrd_genes)) %>%
  filter(!is.element(Gene,brca_genes))%>%
  filter(TVAF > 0.1) %>%
  filter(ExonicFunc!= "synonymous SNV") %>%
  filter(Func == "exonic"|is.element(Func,c("splicing","exonic;splicing")))%>%
  filter(AF < 0.01) %>%
  filter(is.element(ExonicFunc,c("frameshift substitution","stopgain","startloss","."))) -> df.hrd_germline #this list has to be discussed with an expert


ids %>% 
  filter(Material=="wb") %>%
  filter(firstTimepoint_wb == 1)%>%
  mutate(hrd_germline = ifelse(is.element(Patient.ID,df.hrd_germline$Patient.ID),1,0))%>%
  dplyr::select(Patient.ID,hrd_germline) %>% 
  unique -> id.hrd_germline

##Save RData for further use
tempdata <-ls()
rm(list=tempdata[!is.element(tempdata,c("df.hrd_germline","id.hrd_germline","df.brca_germline","id.brca_germline"))])
rm(tempdata)

save.image("data/interim/brca.RData")

#filename="output/BRCA_and_HRD_GermlineStatus.xlsx"
#write.xlsx(df.hrd_germline,filename,sheetName="HRD", append=TRUE)

