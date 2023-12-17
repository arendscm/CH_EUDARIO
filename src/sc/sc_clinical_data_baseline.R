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
source("src/global_functions_themes.R")

########## Load clinical data  ##########
load("data/interim/sc_clin.RData")

########   Create Table with baseline characteristics    ########
#baseline variables
my_vars_baseline=c("Age_atBE",
                   "history.of.cancer",
                   "brca1_germline",
                   "brca2_germline",
                   "brca_germline",
                   "no_prev_lines_binom",
                   "Duration_PriorPARPi",
                   "Duration_PARPi_level")

cat_vars_baseline=c("brca1_germline",
                    "brca2_germline",
                    "brca_germline",
                    "history.of.cancer",
                     "no_prev_lines_binom",
                    "Duration_PARPi_level")

cont_vars_baseline = setdiff(my_vars_baseline,cat_vars_baseline)

df.clin_sc %>% 
  mutate(brca_germline = sign(brca1_germline + brca2_germline)) %>%
           CreateTableOne(strata = "CH",
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
df.clin_sc %>% mutate(nom = ifelse(nom_CH==0,"0",ifelse(nom_CH==1,"1",">1")))%>%
  mutate(nom = factor(nom,levels=c("0","1",">1")))%>%
  ggplot(aes(x = Age_atBE)) + 
  geom_histogram(aes(y=..count..,fill=as.factor(nom)),size=1,position="stack",binwidth=8)  +
  scale_fill_manual(values=scales::seq_gradient_pal(high = "#E64B35FF", low = "#4DBBD5FF",space="Lab")(seq(0,1,length.out=3)), 
                    name="No. of mutations") +
  xlab("Age in years") + 
  ylab("No. of patients") + 
  my_theme()  -> p.agedens 
p.agedens

png("output/figures/sc/sc_agedens.png",width=4, height=3,units="in",res=500,type="cairo")
p.agedens
dev.off()

