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

########   Load IDs    ########
source("src/material_table.R")

########## Load clinical data  ##########
load("data/interim/clin.RData")

########   Create Table with baseline characteristics    ########
#baseline variables
my_vars_baseline=c("Arm",
                   "BRCA1",
                   "BRCA2",
                   "ECOG",
                   "LVEF_C1D1",
                   "TumorBurden_baseline",
                   "Number_PreviousLines",
                   "Number_PreviousPlatinumLines",
                   "No_Platinum_lines_binom",
                   "Type_PreviousTherapy",
                   "PriorPARPi",
                   "Duration_PriorPARPi",
                   "Age_TreatmentStartEUDARIO")

cat_vars_baseline=c("Arm",
                   "BRCA1",
                   "BRCA2",
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
        #nonnormal=cont_vars_baseline,
        exact=cat_vars_baseline,
        missing=TRUE,
        showAllLevels=TRUE,
        quote=FALSE) -> baseline.csv

write.csv(baseline.csv, file = "output/tables/baseline_CH.csv")


