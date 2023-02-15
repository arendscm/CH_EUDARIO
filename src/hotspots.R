# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: List of hotspot mutations described in myeloid malignancies
#
# Input: df from datapreprocessing, list of hotspot variants data/external/hotspots_mm.csv
#
# Output: list of hotspots present in df saved as .RData
#
# ==============================================================================
########   Dependencies   #####
library(base)
library(dplyr)
library(xlsx)
library(reshape)
library(tidyr)
library(readxl)
library(reshape2)
library(stringr)
library(data.table)

######## loading seq data #####
load('data/interim/seqdata.RData')

##alternative hotspot approach: List of hotspots found in myeloid malignancies as described in Feusier et al Blood Cancer Discovery 2021
hotspots <- read.csv("data/external/hotspots_mm.csv",sep=";") %>% 
  mutate(id = paste(Gene,transcript,AA_locus,sep=":"))

##this step takes time!
df%>%
  dplyr::select(AAChange)%>%
  unique%>%
  mutate(AAChange_mod=AAChange)%>%
  separate_rows(AAChange_mod, sep = ",")%>%
  mutate(AAChange_mod = str_replace(AAChange, ":exon.*:p.", ":p."))%>%
  mutate(MM_hotspot = apply(sapply(X = hotspots$id, FUN = grepl, AAChange_mod), MARGIN =  1, FUN = any))%>%
  dplyr::select(AAChange,MM_hotspot)%>% unique %>% filter(MM_hotspot) -> mm_hotspots

##Save RData for further use
  tempdata <-ls()
  rm(list=tempdata[tempdata != "mm_hotspots"])
  rm(tempdata)
  
save.image('data/interim/hotspots.RData')
  