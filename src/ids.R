library(base)
library(dplyr)
library(xlsx)
library(reshape)
library(tidyr)
library(readxl)
library(reshape2)
library(stringr)

##Patient ID table that identifies Sample IDs with Patient ID and timepoints
ids2 <- read_excel("data/external/Ovarial-Ca_LibPrep.xlsx", 
                  sheet = "PatIDfix")
ids2 %>% filter(!is.na(Patient.ID)) -> ids 

#internal Sample IDs
Sample_IDs <- read_excel("data/external/Sample Registry.xlsx", 
                         sheet = "Sample IDs")
select(Sample_IDs,'Sample_orig','Sample.ID', 'External Sample ID', 'Visite', 'Internal Pat ID')->Sample_IDs
left_join(ids,Sample_IDs)->ids

rm(Sample_IDs)
