
# ==============================================================================
# Ovarian Cancer EUDARIO
#
# Author: Max & Klara
#
# Description: Create table with available material at different timepoints for each patient
#
# Input: ids from ids.R
#
# Output: table df.material
# press ALT-O
# ==============================================================================
########   Dependencies   #####
library(base)
library(dplyr)
library(tidyr)
library(fastDummies)
##Patient ID table that identifies Sample IDs with Patient ID and timepoints
source("src/ids.R")

##### create a table with available material at timepoints for each patient ####
ids %>% 
  mutate(visit_mat = paste(Visite,Material,sep="_")) %>% 
  data.frame %>%
  dummy_cols(.,select_columns = "visit_mat") %>% 
  dplyr::select(Patient.ID,
                visit_mat_C1D1_cf,
                visit_mat_C1D1_wb,
                visit_mat_C7D1_cf,
                visit_mat_EOT_cf,
                visit_mat_EOT_wb)%>% 
  group_by(Patient.ID) %>% 
  mutate(c1d1_cf = sum(visit_mat_C1D1_cf),
         c1d1_wb = sum(visit_mat_C1D1_wb),
         c7d1_cf = sum(visit_mat_C7D1_cf),
         eot_cf= sum(visit_mat_EOT_cf),
         eot_wb = sum(visit_mat_EOT_wb))%>%
  data.frame %>%
  dplyr::select(Patient.ID,
                c1d1_cf,
                c1d1_wb,
                c7d1_cf,
                eot_cf,
                eot_wb)%>%
  unique %>%
  mutate(sum_wb = c1d1_wb+eot_wb,
         sum_cf = c1d1_cf + c7d1_cf + eot_cf)-> df.material

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

