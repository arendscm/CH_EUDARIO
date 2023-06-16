# ==============================================================================
# Single Cell Samples Analysis
#
# Author: Max & Klara
#
# Description: Seq Analysis of SC samples and follow up
#
# Input: seqdata
#
# Output: plots
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
library(ggplot2)
library(RColorBrewer)


########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')

load("data/interim/newsamples_filtered.RData")
load('data/interim/seqdata.RData')
SC_registry <- read_excel("data/external/SC_registry.xlsx",sheet = "Patient details")
PatientandSample_Details <- read_excel("data/external/Proben-Registry_Ovarial-CA.xlsx", 
                                       sheet = "Patient details")

################### Overview table ####
#works once tags are included in seqdata preparation
left_join(df.filtered,SC_registry, by="Sample_orig")%>%
  mutate(Patient.ID = Patient.ID.y)->df.filtered
df.filtered%>%
  filter(tag== "true")->Patresults

ch_genes <- c("DNMT3A","TET2","ASXL1","PPM1D","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5")
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")

##Patient.ID
select(df.filtered,Patient.ID)%>%
  unique()->Patient.ID


#filtered results of patients
Patresults%>%
  filter(TVAF >= 0.05)->Patresults5
Patresults%>%
  filter(TVAF >= 0.02)->Patresults2
Patresults%>%
  filter(TVAF >= 0.1)->Patresults10
Patresults%>%
  filter(TVAF >= 0.01)->Patresults1


#Create a new dataframe with patient IDs as the first column
Table <- Patient.ID %>%
  mutate(CH = 0, `VAF >10%` = 0, `VAF >5%` = 0, `VAF >2%` = 0,`VAF >1%` = 0,
         `Mutation_comb+VAF`= NA)


# For each patient in the mutations dataframe, check if they have a mutation in any of the genes
# If they do, add 1 to the corresponding gene column in the new table
for (i in 1:nrow(Patresults)) {
  ID <- Patresults[i, "Patient.ID"]
  ID<-as.character(ID)
  gene <- Patresults[i, "Gene"]
  gene<-as.character(gene)
  if (gene %in% ch_genes) {
    Table[Table$Patient.ID == ID, "CH"] <- 1
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

## Mutation with VAF over 1%
# Get the unique patient ids from the "Patresults2" table
patient_ids_1 <- unique(Patresults1$Patient.ID)

# Compare the patient_ids from Table with the patient_ids in Patresults5
Table$`VAF >1%` <- ifelse(Table$Patient.ID %in% patient_ids_1, 1, Table$`VAF >1%`)

Table[is.na(Table)] <- 0


##Mutation comb with VAF
Patient.ID$Patient.ID%>%
  unique()->ID

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


filename="output/Overview_Table_SC.xlsx"
write.xlsx2(Table, filename,sheetName = "Overview", append=TRUE, overwrite=TRUE)

tempdata <-ls()
rm(list=tempdata[tempdata != "df.filtered"])
rm(tempdata)

##Comutational plots
###Comutational plots
df.filtered%>%
  filter(tag == "true")->test

test<-test%>%
  group_by(Patient.ID) %>%
  filter(n() >= 2) %>%
  ungroup()

Comutations_all <-ggplot(test, aes(x = position, y = TVAF, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Patient.ID, scales = "free_x", ncol = 6) +
  labs(x= "",y = "TVAF", fill = "Genes")+
  theme(axis.text.x = element_blank(),
        xlab = NULL,
        axis.text.y = element_text(size = 8))+
  scale_y_continuous(breaks = c(0.01, 0.1, 0.2, 0.3),
                     labels = c(0.01, 0.1, 0.2, 0.3))
#scale_fill_brewer(palette = "PuBuGn")
Comutations_all

png("output/figures/comutationsSC.png",width=15, height=5,units="in",res=500,type="cairo")
Comutations_all
dev.off()

#### Analysis ####
PatientandSample_Details->data
#Violin Age
ggplot(data, aes(x = ifelse(Multiple_Mutations == 1, 0.2, -0.2), y = Age, fill = factor(Multiple_Mutations), group = Multiple_Mutations)) + 
  geom_violin(alpha = 0.5, position = position_nudge(x = 0.2)) +
  geom_boxplot(aes(x = ifelse(Multiple_Mutations == 1, 0.2, -0.2)), 
               width = 0.05, fill = "white", color = "black", position = position_nudge(x = 0.2)) +  labs(x = "", y = "Age", fill = "Multiple Mutations") +
  scale_fill_manual(values = c("lightblue", "pink"), 
                    labels = c("No", "Yes")) +
  theme_minimal()+
  theme(legend.title = element_text(size = 10))+
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.35, face = "italic", size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))->p.violin.age

png("output/figures/SC-p.violin.age.png",width=10, height=10,units="in",res=500,type="cairo")
p.violin.age
dev.off()

#Violin Months PARPi
ggplot(data, aes(x = ifelse(Multiple_Mutations == 1, 0.2, -0.2), y = months_PARPi, fill = factor(Multiple_Mutations), group = Multiple_Mutations)) + 
  geom_violin(alpha = 0.5, position = position_nudge(x = 0.2)) +
  geom_boxplot(aes(x = ifelse(Multiple_Mutations == 1, 0.2, -0.2)), 
               width = 0.05, fill = "white", color = "black", position = position_nudge(x = 0.2)) +  labs(x = "", y = "Months PARPi Treatment", fill = "Multiple Mutations") +
  scale_fill_manual(values = c("lightblue", "pink"), 
                    labels = c("No", "Yes")) +
  theme_minimal()+
  theme(legend.title = element_text(size = 10))+
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.35, face = "italic", size = 16),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))->p.violin.monthsPARPi

png("output/figures/SC-p.violin.monthsPARPi.png",width=10, height=10,units="in",res=500,type="cairo")
p.violin.monthsPARPi
dev.off()

#Violin Rezidive
ggplot(data, aes(x = ifelse(Multiple_Mutations == 1, 0.2, -0.2), y = Rezidiv, fill = factor(Multiple_Mutations), group = Multiple_Mutations)) + 
  geom_violin(alpha = 0.5, position = position_nudge(x = 0.2)) +
  geom_boxplot(aes(x = ifelse(Multiple_Mutations == 1, 0.2, -0.2)), 
               width = 0.05, fill = "white", color = "black", position = position_nudge(x = 0.2)) +  labs(x = "", y = "Rezidiv", fill = "Multiple Mutations") +
  scale_fill_manual(values = c("lightblue", "pink"), 
                    labels = c("No", "Yes")) +
  theme_minimal()+
  theme(legend.title = element_text(size = 10))+
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.35, face = "italic", size = 16),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))->p.violin.Rezidiv

png("output/figures/SC-p.violin.Rezidiv.png",width=10, height=10,units="in",res=500,type="cairo")
p.violin.Rezidiv
dev.off()

