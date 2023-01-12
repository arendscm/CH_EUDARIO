# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: QC of samples
#
# Input: seqdata, QC fastq files
#
# Output: QC Sequencing runs, Heatmap, Comparison table between germline profiles
#         of all samples one patient has
#
# ==============================================================================
########   Dependencies   #####
library(base)+
library(dplyr)+
library(ggplot2)+
library(xlsx)+
library(stringr)+
library(ggthemes)+
library(viridis)+
library(reshape)+
library(ggpubr)+
library(g3viz)+
library(tidyr)+
library(readxl)+
library(reshape2)

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########   QC sequencing runs  ####
HSM_1346 <- read_excel("QC/QC_runs.xlsx", sheet = "HSMetrics-P1346", 
                      col_types = c("text", "numeric", "numeric"))
GS_1346 <- read_excel("QC/QC_runs.xlsx", sheet = "General Statistics-P1346", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))
HSM_1519 <- read_excel("QC/QC_runs.xlsx", sheet = "HSMetrics-P1519", 
                       col_types = c("text", "numeric", "numeric"))
GS_1519 <- read_excel("QC/QC_runs.xlsx", sheet = "General Statistics-P1346", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))
full_join(HSM_1346,HSM_1519)->HSM
rm(HSM_1346)
rm(HSM_1519)
bind_rows(GS_1346,GS_1519)->GS
rm(GS_1346)
rm(GS_1519)

left_join(HSM,GS,by="Sample Name")->QC
rm(HSM)
rm(GS)
summary(QC)->d
rm(QC)

filename="reports/QC/QC-sequencing runs_means.xlsx"
write.xlsx(d,filename,append=TRUE)
rm(d)

