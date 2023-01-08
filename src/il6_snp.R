# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: determins IL6-SNP status of patients from interim seqdata and ids
#
# Input: data/interim/seqdata.R and src/ids.R 
#
# Output: list of Patient ids with il6snp status
#
# ==============================================================================
########   Dependencies   #####
library(base)
library(dplyr)
library(tidyr)


########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')

######## Get Patient ids
source("src/ids.R")

######## Determin IL6 SNP status of samples
ids %>% filter(!str_detect(Sample_orig,"cf")) %>%
  left_join(.,df %>% 
              filter(Gene == "IL6R") %>% 
              filter(Genomic_Coordinate_hg38 == "chr1:g.154454494:A>C") %>%
              dplyr::select(Sample,Patient.ID,TVAF))%>%
  mutate(il6r_snp = ifelse(TVAF > 0.9,2,
                           ifelse(TVAF < 0.1,0,1))) %>% 
  mutate(il6r_snp = ifelse(is.na(il6r_snp),0,il6r_snp))%>%
  #mutate(Patient.ID = Patient.ID.x)%>%
  dplyr::select(Patient.ID,il6r_snp) %>% 
  unique -> id.il6

rm(df)
rm(ids)

save.image("data/interim/il6snp.RData")