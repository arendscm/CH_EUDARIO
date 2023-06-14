
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
# 
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
  filter(is.na(replicate))%>%
  mutate(visit_mat = paste(Visite,Material,sep="_")) %>% 
  data.frame %>%
  dummy_cols(.,select_columns = "visit_mat") %>% 
  dplyr::select(Patient.ID,
                visit_mat_C1D1_cf,
                visit_mat_C1D1_wb,
                visit_mat_C7D1_cf,
                visit_mat_UE_cf,
                visit_mat_EOT_cf,
                visit_mat_EOT_wb)%>% 
  group_by(Patient.ID) %>% 
  mutate(c1d1_cf = sum(visit_mat_C1D1_cf),
         c1d1_wb = sum(visit_mat_C1D1_wb),
         c7d1_cf = sum(visit_mat_C7D1_cf),
         ue_cf = sum(visit_mat_UE_cf),
         eot_cf= sum(visit_mat_EOT_cf),
         eot_wb = sum(visit_mat_EOT_wb))%>%
  data.frame %>%
  dplyr::select(Patient.ID,
                c1d1_cf,
                c1d1_wb,
                c7d1_cf,
                ue_cf,
                eot_cf,
                eot_wb)%>%
  unique %>%
  mutate(sum_wb = c1d1_wb+eot_wb,
         sum_cf = c1d1_cf + c7d1_cf + ue_cf + eot_cf)-> df.material
