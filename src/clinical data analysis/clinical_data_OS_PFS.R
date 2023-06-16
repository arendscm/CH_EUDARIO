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
library(survival)
library(forestmodel)
library(finalfit)
library(ggplot2)
library(survminer)
library(viridis)

########   Load theme    ########
source("src/global_functions_themes.R")

########## Load clinical data  ##########
load("data/interim/clin.RData")
df -> df.surv

########## Crude OS analysis with Kaplan Meier ##############

surv_obj <- Surv(time = df.surv$OS_days, event = df.surv$OS_event)
fit.km <- survfit(surv_obj ~ CH, data = df.surv)

df.surv %>% 
  ggsurvplot(fit.km, 
             data = ., 
             pval = TRUE,
             #conf.int = TRUE,
             palette = viridis(3),
             risk.table=TRUE,
             tables.height = 0.3,
             ylab = "Overall Survival",
             xlab = "Time in Days",
            # legend.title = "CH status",
            # legend.labs = c("negative", "positive")
  ) -> p.surv
#legend.labs = c("Group1", "Group2"),
# ggtheme = theme_Publication())
p.surv

png("output/figures/surv_CH.png",width=6, height=6,units="in",res=500,type="cairo")
p.surv
dev.off()


##analyse impact of covariates
covariates = c("Age_TreatmentStartEUDARIO","ECOG","PriorPARPi","No_Platinum_lines_binom","TumorBurden_baseline","BRCA1","BRCA2","CH")

finalfit(df.surv, dependent="Surv(OS_days,OS_event)",explanatory=covariates)


######## Crude PFS analysis #################

surv_obj <- Surv(time = df.surv$PFS_days, event = df.surv$PFS_event)
fit.km <- survfit(surv_obj ~ CH, data = df.surv)

df.surv %>% 
  ggsurvplot(fit.km, 
             data = ., 
             pval = TRUE,
             #conf.int = TRUE,
             palette = viridis(3),
             risk.table=TRUE,
             tables.height = 0.3,
             ylab = "Progression-free Survival",
             xlab = "Time in Days",
             #legend.title = "CH status",
             #legend.labs = c("negative", "positive")
  ) -> p.surv
#legend.labs = c("Group1", "Group2"),
# ggtheme = theme_Publication())
p.surv

png("output/figures/surv_CH.png",width=6, height=6,units="in",res=500,type="cairo")
p.surv
dev.off()

covariates = c("Age_TreatmentStartEUDARIO","ECOG","PriorPARPi","No_Platinum_lines_binom","TumorBurden_baseline","BRCA1","BRCA2","CH")

finalfit(df.surv, dependent="Surv(PFS_days,PFS_event)",explanatory=covariates)

