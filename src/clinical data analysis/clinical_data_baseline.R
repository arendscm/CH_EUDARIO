# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: create master tables with clinical data and mutational data
#
# Input: data/extrenal/clindat_modified.xlsx, df.filtered_c1d1
#
# Output: data.frame with clinical and genomic data
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
library(dplyr)
library(tableone)
library(ggplot2)
library(ggsci)

########   Global input    ########
source("src/material_table.R")
source("src/global_functions_themes.R")

########## Load clinical data  ##########
load("data/interim/clin.RData")
df.clin -> df
########   Create Table with baseline characteristics    ########
#baseline variables
my_vars_baseline=c("Age_TreatmentStartEUDARIO",
                   "Arm",
                   "ECOG",
                   "TumorBurden_baseline",
                   "LVEF_C1D1",
                   "BRCA1",
                   "brca1_germline",
                   "BRCA2",
                   "brca2_germline",
                   "Number_PreviousLines",
                   "Number_PreviousPlatinumLines",
                   "No_Platinum_lines_binom",
                   "Type_PreviousTherapy",
                   "PriorPARPi",
                   "Duration_PriorPARPi")

cat_vars_baseline=c("Arm",
                    "BRCA1",
                    "brca1_germline",
                    "BRCA2",
                    "brca2_germline",
                   "ECOG",
                   "Type_PreviousTherapy",
                   "PriorPARPi",
                   "No_Platinum_lines_binom")

cont_vars_baseline = setdiff(my_vars_baseline,cat_vars_baseline)

df %>% CreateTableOne(strata = "CH",
                      vars=c(my_vars_baseline),
                      factorVars = cat_vars_baseline,
                      includeNA=FALSE,
                      #addOverall = TRUE,
                      data=.) %>% 
  print(., 
        nonnormal=cont_vars_baseline,
        exact=cat_vars_baseline,
        missing=TRUE,
        showAllLevels=TRUE,
        quote=FALSE) -> baseline.csv

write.csv(baseline.csv, file = "output/tables/baseline_CH.csv")
write.xlsx(baseline.csv, file = "output/tables/baseline.xlsx",sheetName = "baseline")


######### Plot age distribution #######################
df %>% mutate(nom = ifelse(nom==0,"0",ifelse(nom==1,"1",">1")))%>%
  mutate(nom = factor(nom,levels=c("0","1",">1")))%>%
  ggplot(aes(x = Age_TreatmentStartEUDARIO)) + 
  geom_histogram(aes(y=..count..,fill=as.factor(nom)),size=1,position="stack",binwidth=10)  +
  scale_fill_npg(breaks=c("0","1",">1"),name="No. mutations") +
  xlab("Age in Years") + 
  #xlim(15,85) +
  ylab("Number of Patients") + 
  #ggtitle("Age Distribution according to CH status") +
  my_theme()  -> p.agedens 
p.agedens

png("output/figures/agedens.png",width=6, height=5,units="in",res=500,type="cairo")
p.agedens
dev.off()

