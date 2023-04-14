# ______________________________________________________________________________
# Ovarian Cancer EUDARIO
#
# Author: Max & Klara
#
# Description: Sample registry and availability
#
# Input: sample registry sheet
#
# Output: sample registry and sample plot
# 
# ______________________________________________________________________________
########   Dependencies   #####
library(base)
library(dplyr)
library(tidyr)
library(readxl)
library(xlsx)

##### Sample Availability  ####
Sample_registry <- read_excel("data/external/Sample Registry.xlsx", sheet = "Sample_registry_Overview", col_types = c("text", "numeric", "text", "text"))

Sample_registry%>%
  filter(Material== "wb" | Material == "cf")->Sample_registry_sq
# Calculate sample counts
sample_counts <- Sample_registry_sq %>% group_by(Timepoint, Material) %>% summarize(count = n())

# Create plot all Samples
ggplot(Sample_registry_sq, aes(x = Timepoint, y = as.numeric(factor(Patient)), color = Material)) + 
  geom_point(position = position_dodge(width = 0.5), size = 1) + 
  scale_color_manual(values = c("blue","red")) + 
  labs(x = "Timepoints", y = "Patients") +
  scale_y_continuous(breaks = seq_along(unique(Sample_registry$Patient)), labels = seq_along(unique(Sample_registry$Patient)))+
  theme(axis.text.y = element_text(size = 5))+
  ggtitle("Samples Sequenced") ->p.Samples_sq

png("output/figures/p.Samples_sq.png",width=7, height=8,units="in",res=500,type="cairo")
p.Samples_sq
dev.off()

##all Samples
# Calculate sample counts
sample_counts <- Sample_registry %>% group_by(Timepoint, Material) %>% summarize(count = n())

# Create plot all Samples
ggplot(Sample_registry, aes(x = Timepoint, y = as.numeric(factor(Patient)), color = Material)) + 
  geom_point(position = position_dodge(width = 0.5), size = 1) + 
  scale_color_manual(values = c("blue", "orange", "grey", "red")) + 
  labs(x = "Timepoints", y = "Patients") +
  scale_y_continuous(breaks = seq_along(unique(Sample_registry$Patient)), labels = seq_along(unique(Sample_registry$Patient)))+
  theme(axis.text.y = element_text(size = 5))+
  ggtitle("All Samples") ->p.Samples

png("output/figures/p.Samples.png",width=7, height=8,units="in",res=500,type="cairo")
p.Samples
dev.off()



##### Sample registry ####
#Sample Registry -> what timepoint/sample do we have from each patient
Ovarial_Ca_LibPrep <- read_excel("data/external/Ovarial-Ca_LibPrep.xlsx", 
                                 sheet = "Pooling", skip=1)
source("src/ids.R")
select(Ovarial_Ca_LibPrep,'External Pat ID', 'Visite','Internal Sample ID', 'External Sample ID')->pool

select(ids,Sample.ID,'External Sample ID', replicate, Int.Patient.ID)->id
full_join(id,pool,by="External Sample ID")%>%
  filter(is.na(replicate))->pool

pool%>%
  filter(Visite=="C1D1-wb")->C1D1
pool%>%
  filter(Visite=="EOT-wb")->EOT
pool%>%
  filter(Visite=="C1D1-cf")->cfC1D1
pool%>%
  filter(Visite=="C7D1-cf")->cfC7D1
pool%>%
  filter(Visite=="EOT-cf")->cfEOT
pool%>%
  filter(Visite=="UE-cf")->cfUE
pool%>%
  filter(Visite == "CxD1-cf")->cfCxD1

full_join(C1D1,EOT, by='External Pat ID') %>%
  full_join(., cfC1D1, by= 'External Pat ID')%>%
  full_join(., cfC7D1, by= 'External Pat ID')%>%
  full_join(., cfEOT, by= 'External Pat ID')%>%
  full_join(.,cfUE,by='External Pat ID')%>%
  full_join(.,cfCxD1, by='External Pat ID')->list
list%>%
  unique()->list

filename="data/external/Sample Registry_test.xlsx"
write.xlsx(list,filename, sheetName="New Sample registry", append=TRUE)
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
                     sheet = "Int-Ext ID")
komb%>%
  mutate(Patient.ID = as.character(Patient.ID))->komb
IntExt%>%
  mutate(Patient.ID = as.character(Patient.ID))->IntExt
left_join(komb,IntExt, by='Patient.ID')->komb2

filename="test3.xlsx"
write.xlsx(komb2, filename, sheetName = "list", append=TRUE)
rm(komb2)
rm(IntExt)
rm(komb)
rm(filename)


