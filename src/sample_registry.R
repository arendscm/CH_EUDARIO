# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: creates table of all samples (Patient, timepoints...) 
#              creates tables with sample IDs ( external and internal ones)
#
# Input: external data folder
#
# Output: table of samples and sample IDs
#
# ==============================================================================
########   Dependencies   #####
library (readxl)
library (dplyr)

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########   Sample Registry and Int-Ext Pat ID ####

#Sample Registry -> what timepoint/sample do we have from each patient
Ovarial_Ca_LibPrep <- read_excel("data/external/Ovarial-Ca_LibPrep.xlsx", 
                                 sheet = "Pooling", skip=1)
select(Ovarial_Ca_LibPrep,'External Pat ID', 'Visite','Internal Pat ID', 'External Sample ID')->pool

pool%>%
  filter(Visite=="C1D1")->C1D1
pool%>%
  filter(Visite=="EOT")->EOT
pool%>%
  filter(Visite=="cf-C1D1")->cfC1D1
pool%>%
  filter(Visite=="cf-C7D1")->cfC7D1
pool%>%
  filter(Visite=="cf-EOT")->cfEOT

full_join(C1D1,EOT, by='External Pat ID') %>%
full_join(., cfC1D1, by= 'External Pat ID')%>%
full_join(., cfC7D1, by= 'External Pat ID')%>%
full_join(., cfEOT, by= 'External Pat ID')->list

filename="data/interim/sample_registry.csv"
write.csv(list, filename)
rm(C1D1)
rm(cfC1D1)
rm(cfC7D1)
rm(cfEOT)
rm(EOT)
rm(pool)
rm(list)
rm(Ovarial_Ca_LibPrep)


#combine Ext and Int PatID
#komb aus Extraction Plan 
komb <- read_excel("data/external/Sample Registry.xlsx", 
                   sheet = "Kombinieren", col_types = c("text", 
                                                        "text", "text", "text"))
#IntExt aus sample registry
IntExt <- read_excel("data/external/Sample Registry.xlsx", 
                     sheet = "Sample IDs")
left_join(komb,IntExt, by='External Pat ID')->komb2

filename="Sample Registry/IntExtPatID.xlsx"
#write.xlsx(komb2, filename, sheetName = "list", append=TRUE)
rm(komb2)
rm(IntExt)
rm(komb)
rm(filename)

