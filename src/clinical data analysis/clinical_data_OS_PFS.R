# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: create master tables with clinical data and mutational data
#
# Input: preprocessed clinical data and mutational data
#
# Output: tables plots for response, PFS and OS
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
library(survival)
library(forestmodel)
library(finalfit)
library(ggplot2)
library(survminer)
library(ggsci)

########   Load theme    ########
source("src/global_functions_themes.R")

########## Load clinical data  ##########
load("data/interim/clin.RData")
df.clin %>% mutate(brca_germline=brca1_germline+brca2_germline) %>%
  mutate(age_median = ifelse(Age_TreatmentStartEUDARIO >= 62,">=62","<62"))-> df.surv

########## Response assessment ######################
my_vars_response=c("Response_best","response_binom","response_binom2")

cat_vars_response=my_vars_response

df.clin %>% CreateTableOne(strata = "CH",
                      vars=c(my_vars_response),
                      factorVars = cat_vars_response,
                      includeNA=FALSE,
                      #addOverall = TRUE,
                      data=.) %>% 
  print(., 
        #nonnormal=cont_vars_ae,
        exact=cat_vars_response,
        missing=TRUE,
        showAllLevels=TRUE,
        quote=FALSE) 

########## Crude OS analysis with Kaplan Meier ##############

surv_obj <- Surv(time = df.surv$OS_days, event = df.surv$OS_event)
fit.km <- survfit(surv_obj ~ CH, data = df.surv)

df.surv %>% 
  ggsurvplot(fit.km, 
             data = ., 
             pval = TRUE,
             #conf.int = TRUE,
             palette = pal_npg("nrc")(2),
             risk.table=TRUE,
             tables.height = 0.3,
             ylab = "Overall Survival",
             xlab = "Time in Days",
            # legend.title = "CH status",
            # legend.labs = c("negative", "positive")
  ) -> p.os

p.os

png("output/figures/surv_CH.png",width=6, height=6,units="in",res=500,type="cairo")
p.os
dev.off()

##analyse impact of covariates
covariates = c("Age_TreatmentStartEUDARIO","ECOG_binom","PriorPARPi","No_Platinum_lines_binom","TumorBurden_baseline","brca_germline","CH")

finalfit(df.surv, dependent="Surv(OS_days,OS_event)",explanatory=covariates)

fit.coxph <- coxph(Surv(time = df.surv$OS_days, event = df.surv$OS_event) ~ CH + Arm + age_median+ TumorBurden_baseline+ ECOG_binom + Number_PreviousLines, 
                   data = df.surv)
summary(fit.coxph)
forest_model(fit.coxph)

######## Crude PFS analysis #################

surv_obj <- Surv(time = df.surv$PFS_days, event = df.surv$PFS_event)
fit.km <- survfit(surv_obj ~ CH, data = df.surv)

df.surv %>% 
  ggsurvplot(fit.km, 
             data = ., 
             pval = TRUE,
             #conf.int = TRUE,
             palette = pal_npg("nrc")(2),
             risk.table=TRUE,
             tables.height = 0.3,
             ylab = "Progression-free Survival",
             xlab = "Time in Days",
             #legend.title = "CH status",
             #legend.labs = c("negative", "positive")
  ) -> p.pfs

p.pfs

png("output/figures/surv_CH.png",width=6, height=6,units="in",res=500,type="cairo")
p.pfs
dev.off()

##analyse impact of covariates
covariates = c("Age_TreatmentStartEUDARIO","ECOG_binom","PriorPARPi","No_Platinum_lines_binom","TumorBurden_baseline","brca1_germline","brca2_germline","CH")

finalfit(df.surv, dependent="Surv(PFS_days,PFS_event)",explanatory=covariates)
#set of covariates that are predictive at p<0.1
covariates_01 = c("ECOG_binom","PriorPARPi","TumorBurden_baseline","CH")
finalfit(df.surv, dependent="Surv(PFS_days,PFS_event)",explanatory=covariates_01)


fit.coxph <- coxph(Surv(time = df.surv$PFS_days, event = df.surv$PFS_event) ~ CH + Arm + age_median+ TumorBurden_baseline+ ECOG_binom + Number_PreviousLines, 
                   data = df.surv)
summary(fit.coxph)
forest_model(fit.coxph)
