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
library(base)
library(dplyr)
library(xlsx)
library(readxl)


########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########   QC sequencing runs  ####
HSM_1346 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "HSMetrics-P1346", 
                      col_types = c("text", "numeric", "numeric"))
GS_1346 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "General Statistics-P1346", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))
HSM_1519 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "HSMetrics-P1519", 
                       col_types = c("text", "numeric", "numeric"))
GS_1519 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "General Statistics-P1519", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))
HSM_1803 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "HSMetrics-P1803", 
                       col_types = c("text", "numeric", "numeric"))
GS_1803 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "General Statistics-P1803", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))


left_join(HSM_1346,GS_1346)->QC_1346
left_join(HSM_1519,GS_1519)->QC_1519
left_join(HSM_1803,GS_1803)->QC_1803


for (file in c("QC_1803","QC_1519","QC_1346")){
  print(file)
  eval(as.name(file))->QC
  #remove samples that didnt work-> "OvCA_45_cf_C1D1.realigned" "OvCA_46_cf_C1D1.realigned" "OvCA_54_cf_C1D1.realigned"
  #QC%>%
    #filter(`Target Bases 30X` <= 0.95)->failed
  #failed$`Sample Name`->failed
  #QC%>%
    #filter(!is.element(`Sample Name`, failed))->QC
  summary(QC)->d
  rm(QC)
  
  filename="output/qc/QC-sequencing runs_means.xlsx"
  write.xlsx(d,filename,sheetName=file,append=TRUE)
  rm(d)
}


