# ______________________________________________________________________________
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: create workingtable->all wb true with corresponding timepoints and Materials
#                                 ->all cf onyl with corresponing timepoints
#
# Input: df from data/interim/seqdata.RData
#
# Output: working table
# ______________________________________________________________________________
#####  Dependencies   #####
library(dplyr)
library(xlsx)
load('data/interim/seqdata_filtered.RData')
load('data/interim/seqdata.RData')
load('data/interim/seqdata_filtered_cf.RData')

#WB
df.filtered.c1d1%>%
  filter(tag != "false")%>%
  filter(tag != "germline")%>%
  .$Patmut->a
df%>%
  filter(is.element(Patmut,a))->workingdata1

#cf
df%>%
  filter(Material == "cf")%>%
  filter(tag != "false")%>%
  filter(tag != "germline")%>%
  filter(!is.na(tag))%>%
  .$Patmut->b
df%>%
  filter(is.element(Patmut,b))->workingdata2
bind_rows(workingdata1, workingdata2)%>%
  filter(is.na(replicate))%>%
  unique()->workingdata

filename <- paste("output/workingdata",Sys.Date(),".xlsx",sep="")
write.xlsx(workingdata, filename, sheetName="all", append=TRUE)


tempdata <-ls()
rm(list=tempdata[tempdata != "workingdata"])
rm(tempdata)

