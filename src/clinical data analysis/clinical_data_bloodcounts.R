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
library(ggplot2)
library(ggsci)


########   Load IDs    ########
source("src/material_table.R")
source("src/genegroup_definitions.R")

########## Load mutation data with clinical data ##########
load("data/interim/seqdata_filtered.RData")

########   Load clinical data    ########
df.bc <- data.frame(read.table("data/external/bc_modified.csv",header=TRUE,sep=";"))


load("data/interim/clin.RData")


### plot blood counts during follow up


df.bc%>%
  filter(phase=="M",day==1)%>%
  ggplot(aes(y=wbc,x=cycle,group=patient_code)) + 
  geom_point(size=1,na.rm=FALSE) + 
  geom_line(size=0.5,na.rm=FALSE) + 
  facet_wrap(~ patient_code, ncol=18, dir="h") +
  scale_y_log10()+
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  scale_color_npg()+
  theme_minimal()-> p.serial_ch
p.serial_ch



