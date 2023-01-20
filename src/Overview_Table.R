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
  filter(Patient.ID != 0)%>%
  unique()->Patient.ID

#filtered results of patients
df.filtered%>%
  filter(Visite == "C1D1")%>%
  filter(TVAF >= 0.005)%>%
  filter(cf == "0")%>%
  filter(tag == "true")%>%
  filter(Patient.ID != 0)->Patresults
Patresults%>%
  filter(TVAF >= 0.05)->Patresults5
Patresults%>%
  filter(TVAF >= 0.02)->Patresults2
Patresults%>%
  filter(TVAF >= 0.1)->Patresults10

#Create a new dataframe with patient IDs as the first column
Table <- Patient.ID %>%
  rename(Patient.ID == "Patient.ID") %>%
  mutate(CH = 0, `VAF >10%` = 0, `VAF >5%` = 0, `VAF >2%` = 0, `only cf-TP53`=0,
         `Mutation_comb+VAF`= NA)
##CH genes ->list of ch_genes
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
##HRD genes ->list of hrd_genes
# Add columns for each gene, with initial values of 0
for (gene in hrd_genes) {
  Table[gene] <- 0
}

# For each patient in the mutations dataframe, check if they have a mutation in any of the genes
# If they do, add 1 to the corresponding gene column in the new table
for (i in 1:nrow(Patresults)) {
  Patient.ID <- Patresults[i, "Patient.ID"]
  gene <- Patresults[i, "Gene"]
  if (gene %in% hrd_genes) {
    Table[Table$Patient.ID == Patient.ID, gene] <- Table[Table$Patient.ID == Patient.ID, gene] + 1
    Table[Table$Patient.ID == Patient.ID, "HRD"] <- 1
  }
}

## Mutation with VAF over 10%
# Get the unique patient ids from the "Patresults5" table
patient_ids_10 <- unique(Patresults10$Patient.ID)

# Compare the patient_ids from Table with the patient_ids in Patresults5
Table$`VAF >10%` <- ifelse(Table$Patient.ID %in% patient_ids_10, 1, Table$`VAF >10%`)

## Mutation with VAF over 5%
# Get the unique patient ids from the "Patresults5" table
patient_ids_5 <- unique(Patresults5$Patient.ID)

# Compare the patient_ids from Table with the patient_ids in Patresults5
Table$`VAF >5%` <- ifelse(Table$Patient.ID %in% patient_ids_5, 1, Table$`VAF >5%`)

## Mutation with VAF over 2%
# Get the unique patient ids from the "Patresults2" table
patient_ids_2 <- unique(Patresults2$Patient.ID)

# Compare the patient_ids from Table with the patient_ids in Patresults5
Table$`VAF >2%` <- ifelse(Table$Patient.ID %in% patient_ids_2, 1, Table$`VAF >2%`)



Table[is.na(Table)] <- 0

genes_total <- c(ch_genes,hrd_genes)

##Mutation comb with VAF
for (patient_id in ID) {
  # only Patient ID results
  patient_mutations <- subset(Patresults, Patresults$Patient.ID == patient_id)
  
  patient_mutations_vaf <- character(nrow(patient_mutations))
  for (i in 1:nrow(patient_mutations)){
    patient_mutations_vaf[i] <- paste0(patient_mutations$Gene[i], " (",patient_mutations$TVAF[i],")")
  }
  # Check if the patient has a mutation
  if (nrow(patient_mutations) > 0) {
    # If the patient has a mutation -> mutation names and VAF
    Table$`Mutation_comb+VAF`[Table$Patient.ID == patient_id] <- paste(patient_mutations_vaf, collapse = " + ")
  }
}

#il6snp
left_join(Table,id.il6, by="Patient.ID")->Table
#BRCA_germine
left_join(Table,id.brca_germline, by="Patient.ID")->Table

#cf Tp53
load('data/interim/df.cf_only.RDATA')
df.cf_only%>%
  filter(Gene == "TP53")%>%
  filter(Sample != "cf2-A7")%>%
  filter(TVAF >=0.01 & TVAF <= 0.5)->cf_TP53

##TP53 Mutation in cf ONLY! (VAF 0.01 - 0.5)
# Get the unique patient ids from the "cf_TP53" table
patient_ids_cf_TP53 <- unique(cf_TP53$Patient.ID)

# Compare the patient_ids from Table with the patient_ids in cf_TP53
Table$`only cf-TP53` <- ifelse(Table$Patient.ID %in% patient_ids_cf_TP53, 1, Table$`only cf-TP53`)


#ids
ids%>%
  filter(Visite == "C1D1")%>%
  select(.,Patient.ID,`Internal Pat ID`)->ids_C1D1
left_join(Table, ids_C1D1,by="Patient.ID")->Table



filename="output/Overview_Table.xlsx"
file.remove("output/Overview_Table.xlsx")
write.xlsx2(Table, filename,sheetName = "Overview", append=TRUE, overwrite=TRUE)
save(Table,file="data/interim/Overview_Table.RDATA")

rm(list = ls())

