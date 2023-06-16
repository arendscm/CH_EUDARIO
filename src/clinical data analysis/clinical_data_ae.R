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
my_vars_ae=c("ae_haematotox",
             "ae_neutropenia",
             "ae_thrombocytopenia",
             "ae_anemia",
             "ae_haematotox_severe")

cat_vars_ae=my_vars_ae

df %>% CreateTableOne(strata = "CH",
                      vars=c(my_vars_ae),
                      factorVars = cat_vars_ae,
                      includeNA=FALSE,
                      #addOverall = TRUE,
                      data=.) %>% 
  print(., 
        #nonnormal=cont_vars_ae,
        exact=cat_vars_ae,
        missing=TRUE,
        showAllLevels=TRUE,
        quote=FALSE) -> ae.csv

write.csv(ae.csv, file = "output/tables/ae_CH.csv")