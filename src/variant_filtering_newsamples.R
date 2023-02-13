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
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')

###include only new samples = those with Patient.ID = NA
df -> df.orig
df.orig %>% filter(is.na(Patient.ID)) -> df

########   FILTERING CH calls------------------------------------------------------------
# filter criteria
#functional criteria
df %>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  dplyr::select(mutID) %>% mutate(mutID = as.character(mutID)) -> mutID.func

## read count criteria
df %>%
  filter(readDepth > 50, TR2 > 9, TVAF > 0.01) %>%
  dplyr::select(mutID) -> mutID.count

## quality criteria
df %>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  dplyr::select(mutID)-> mutID.qual

## mutation call frequency criteria
df %>% 
  filter(AF<0.1) %>%    #filtert alle h?ufig in Datenbanken (=seq errors) gelisteten mutation raus
  filter((mutFreq < 10)|((p.binom<= -10)&med.vaf < 0.44))%>% #filtert alle mutationen raus, die in mehr als 20% der samples auftreten
  dplyr::select(mutID) -> mutID.freq

# rescue ASXL1 dupG mutations <- this step is no longer needed, when we use p.binom 
#df %>%
#  filter(str_detect(Start,'32434638'))%>%filter(str_detect(AAChange,'ASXL1')) %>%
#  mutate(dev.med = ((TVAF - median(TVAF))/sd(TVAF))) %>% #calculate deviation from median in terms of standarddeviations
#  filter(dev.med > 1) %>%    #to be discussed
#  dplyr::select(mutID) -> mutID.asxl1

## rescue hotspots
df %>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 3) %>%
  filter(TVAF >0.005) %>%
  filter(hotspot)%>%
  dplyr::select(mutID)-> mutID.hotspots

## rescue known CHIP mutations with weakened quality criteria
df %>% filter(ChipPub != "") -> mutID.CHIP
df %>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 7) %>%
  filter(TVAF >0.005) %>%
  dplyr::select(mutID)-> mutID.CHIP.qual


#filtering
##somatic variants
inner_join(mutID.func,mutID.count) %>% 
  inner_join(.,mutID.freq) %>% 
  inner_join(.,mutID.qual) %>% 
  inner_join(.,df) %>% 
  full_join(.,inner_join(df,inner_join(mutID.CHIP,mutID.CHIP.qual))) %>%
  #full_join(.,inner_join(df,mutID.asxl1)) %>%
  full_join(.,inner_join(df,mutID.hotspots))%>%
  filter(snp == FALSE) %>%
  mutate(current_filter = 1) %>% ##tag all variants passing current filter, then join with list of previously tagged true
  full_join(df %>% filter(tag=="true"))-> df.filtered


rm(mutID.CHIP)+
rm(mutID.CHIP.qual)+
rm(mutID.count)+
rm(mutID.freq)+
rm(mutID.func)+
rm(mutID.hotspots)+
rm(mutID.qual)+
rm(ids)+
rm(df)

save.image("data/interim/newsamples_filtered.RData")

filename <- paste("output/filtered_newsample_",Sys.Date(),".xlsx",sep="")
write.xlsx(df.filtered,filename,sheetName = "filtered_results",append=TRUE)

