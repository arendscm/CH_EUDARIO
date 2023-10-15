# ==============================================================================
# Single Cell Seq Data Analysis
#
# Author: Max & Klara
#
# Description: Analysis of single cell seq data
#
# Input: XX
#
# Output: plots
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
library(ggplot2)
library(ggsci)


########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')

load("data/interim/newsamples_filtered.RData")
load('data/interim/seqdata.RData')
SC_registry <- read.xlsx("data/external/SC_registry.xlsx",sheetIndex=1,header=TRUE)

######## Functions and themes
source("src/global_functions_themes.R")
source("src/genegroup_definitions.R")

####### load single cell seq data
##filenames
filenames <- list.files(path = "data/raw/ClusterResults/results", pattern = "gt.csv" , all.files = FALSE,
                        full.names = FALSE, recursive = FALSE)

##load files into list of data.frames
all_csv <- lapply(paste("data/raw/ClusterResults/results/",filenames,sep=""),read.csv)

#extract sample names
sample_names <- data.frame(filenames) %>% 
  mutate(file=filenames)%>%
  separate(filenames, c("Sample1","Sample2","dummy","whichSample","gt"),sep="_") %>% 
  mutate(Sample=ifelse(whichSample==0,Sample1,Sample2))%>%
  dplyr::select(-dummy,-gt,-Sample1,-Sample2,-whichSample) 

#names list entries by sample name
names(all_csv) <- sample_names$Sample

#genotype frequency
all_csv[[10]] %>% dplyr::select(-cluster_id,-X) %>% unite(col = "gt",sep="_",remove=FALSE) %>% dplyr::select(gt) %>% table 
