# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: creates overview table of all patients
#
# Input: seqdata, filtered data, BRCA Status, il6snp status
#
# Output: Overview table
# ==============================================================================
########   Dependencies   #####
library(dplyr)
library(tidyr)
library(xlsx)



########  Load preprocessed sequencing data, filtered results,...
load('data/interim/seqdata.RData')
load('data/interim/seqdata_filtered.RData')
load('data/interim/brca.RData')
load('data/interim/il6snp.RData')

######## Get Patient ids
source("src/ids.R")

##Patient.ID
select(df,Patient.ID)%>%
  unique()->Patient.ID

#filtered results of patients
df.filtered%>%
  filter(Visite == "C1D1")%>%
  filter(TVAF >= 0.005)%>%
  filter(cf == "0")%>%
  filter(tag == "true")->Patresults

#Create a new dataframe with patient IDs as the first column
Table <- Patient.ID %>%
  rename(Patient.ID == "Patient.ID") %>%
  mutate(CH = 0, il6snp = 0)

##CH genes
# Add columns for each gene, with initial values of 0
  for (gene in ch_genes) {
    Table[gene] <- 0
  }
  
# For each patient in the mutations dataframe, check if they have a mutation in any of the genes
# If they do, add 1 to the corresponding gene column in the new table
  for (i in 1:nrow(Patresults)) {
    Patient.ID <- Patresults[i, "Patient.ID"]
    gene <- Patresults[i, "Gene"]
    if (gene %in% ch_genes) {
      Table[Table$Patient.ID == Patient.ID, gene] <- Table[Table$Patient.ID == Patient.ID, gene] + 1
      Table[Table$Patient.ID == Patient.ID, "CH"] <- 1
    }
  }

Table[is.na(Table)] <- 0

# Create a new column called "Mutation_comb"
Table$Mutation_comb <- ""

# For each patient, find the genes they have mutations in and concatenate them together with "+"
for (i in 1:nrow(Table)) {
  patient_genes <- c()
  for (gene in ch_genes) {
    if (Table[i, gene] > 0) {
      patient_genes <- c(patient_genes, gene)
    }
  }
  Table[i, "Mutation_comb"] <- str_c(patient_genes, collapse = "+")
}

#il6snp
left_join(Table,id.il6, by="Patient.ID")->Table
#BRCA_germine
left_join(Table,id.brca_germline, by="Patient.ID")->Table

#ids
ids%>%
  filter(Visite == "C1D1")%>%
  select(.,Patient.ID,`Internal Pat ID`)->ids_C1D1
left_join(Table, ids_C1D1,by="Patient.ID")->Table



filename="output/Overview_Table.xlsx"
write.xlsx(Table, filename, sheetName="Overview", append=TRUE)
save(Table,file="data/interim/Overview_Table.RDATA")

rm(list = ls())

