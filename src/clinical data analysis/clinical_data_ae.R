# ==============================================================================
# CH in EUDARIO
#
# Author: Max & Klara
#
# Description: Assess occurence of adverse events during study for CH+ and CH- patients
#
# Input: preprocessed clinical data
#
# Output: table with AEs 
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
             "ae_haematotox_severe",
             "ae_infection",
             "ae_infection_severe",
             "ae_allergic_reaction",
             "ae_allergic_reaction_severe",
             "ae_bleeding",
             "ae_kidney_failure",
             "ae_transaminases")

my_vars_ae_ctx=c("ae_haematotox_ctx",
             "ae_neutropenia_ctx",
             "ae_thrombocytopenia_ctx",
             "ae_anemia_ctx",
             "ae_haematotox_severe_ctx",
             "ae_infection_ctx",
             "ae_infection_severe_ctx",
             "ae_allergic_reaction_ctx",
             "ae_allergic_reaction_severe_ctx",
             "ae_bleeding_ctx",
             "ae_kidney_failure_ctx",
             "ae_transaminases_ctx")

my_vars_ae_main=c("ae_haematotox_main",
                 "ae_neutropenia_main",
                 "ae_thrombocytopenia_main",
                 "ae_anemia_main",
                 "ae_haematotox_severe_main",
                 "ae_infection_main",
                 "ae_infection_severe_main",
                 "ae_allergic_reaction_main",
                 "ae_allergic_reaction_severe_main",
                 "ae_bleeding_main",
                 "ae_kidney_failure_main",
                 "ae_transaminases_main")

##overall AEs
cat_vars_ae=my_vars_ae

df.clin %>% 
  CreateTableOne(strata = "CH",
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

##AEs during Ctx

cat_vars_ae=my_vars_ae_ctx

df.clin %>% 
  CreateTableOne(strata = "DNMT3A",
                 vars=c(my_vars_ae_ctx),
                 factorVars = cat_vars_ae,
                 includeNA=FALSE,
                 #addOverall = TRUE,
                 data=.) %>% 
  print(., 
        #nonnormal=cont_vars_ae,
        exact=cat_vars_ae,
        missing=TRUE,
        showAllLevels=TRUE,
        quote=FALSE) 


##AEs during maintenance
cat_vars_ae=my_vars_ae_main

df.clin %>% 
  filter(Date_StartMaintenanceEUDARIO != "NAP")%>%
  CreateTableOne(strata = "CH",
                 vars=c(my_vars_ae_main),
                 factorVars = cat_vars_ae,
                 includeNA=FALSE,
                 #addOverall = TRUE,
                 data=.) %>% 
  print(., 
        #nonnormal=cont_vars_ae,
        exact=cat_vars_ae,
        missing=TRUE,
        showAllLevels=TRUE,
        quote=FALSE) 


log.reg <- glm(ae_infection ~  DDR + CH, family="binomial",data=df.clin %>% mutate(age_dec=Age_TreatmentStartEUDARIO/10))
summary(log.reg)
