# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: Analysis of longitudinal data
#
# Input: df from data/interim/seqdata.RData
#
# Output: Excel list of filtered results, plots, ...
#
# ==============================================================================
########   Dependencies   #####
library(base)
library(dplyr)
library(stringr)
library(reshape)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(viridis)
library(reshape)
library(ggpubr)
library(g3viz)
library(maftools)

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')
load('data/interim/seqdata_filtered.RDATA')

######## Get Patient ids
source("src/ids.R")

######## Functions and themes

source("src/global_functions_themes.R")

########   SERIAL SAMPLES ####

#find Clones found in C1D1 and EOT for each patient
df.filtered%>%
  filter(Visite== "EOT")%>%
  select(.,Patient.ID)%>%
  filter(Patient.ID != "0")%>%
  unique()%>%
  unlist()->EOT_ids


#für ersten Durchlauf -> table zur Verfügung stellen
Dynamicallpat<-Dynamic

for(current_patient in EOT_ids)
{print(current_patient)
  #Datensätze auf Pat anpassen 
  df%>%
    filter(Patient.ID==current_patient)%>%
    filter(mutFreq<=15)%>%
    filter(snp=="FALSE")%>%
    filter(Visite != "cf-C1D1")->df.pat
  df.pat%>%
    filter(Visite=="C1D1")->df.c1d1pat
  df.pat%>%
    filter(Visite=="EOT")->df.eotpat
  
  #Überschneidungen finden
  inner_join(df.c1d1pat,df.eotpat,by="Patmut")->C1D1EOTdynamic
  
  #Matching Mutationen markieren 
  C1D1EOTdynamic%>%
    mutate(match="yes")->C1D1EOTdynamic
  select(C1D1EOTdynamic,Patmut,match)->Dynamicmut
  
  #mutationen aus dem großen pool patpool rausholen 
  left_join(df.pat,Dynamicmut,by="Patmut")->Dynamic
  Dynamic%>%
    filter(match=="yes")->Dynamic
  
  ##bis jetzt sind einfach alle matching mutationen mit angewandten Filtern drin, viele sind aber sehr klein!
  ##jetzt: C1D1 >0,5 filtern und dann wieder im EOT pool gucken, ob es matches gibt!
  Dynamic%>%
    filter(Visite=="C1D1")%>%
    filter(TVAF >= 0.005)->DynamicC
  Dynamic%>%
    filter(Visite=="EOT")->DynamicE
  
  #Überschneidungen finden
  inner_join(DynamicC,DynamicE,by="Patmut")->C1D1EOTdynamic
  
  #Matching Mutationen markieren 
  C1D1EOTdynamic%>%
    mutate(match2="yes")->C1D1EOTdynamic
  select(C1D1EOTdynamic,Patmut,match2)->Dynamicmut
  
  #mutationen aus dem großen pool patpool rausholen 
  left_join(Dynamic,Dynamicmut,by="Patmut")->Dynamicx
  Dynamicx%>%
    filter(match2=="yes")->Dynamicx
  
  ##Jetzt müssen noch EOT >0.5 gerettet werden die in C1D1 noch kein match hatten!
  DynamicE%>%
    filter(TVAF >=0.005)->DynamicEresc
  #die, die noch kein match haben
  anti_join(DynamicEresc,DynamicC, by="Patmut")->DynamicEresc
  #
  DynamicEresc%>%
    mutate(match3="yes")->DynamicEresc
  select(DynamicEresc,Patmut,match3)->Dynamicmut
  # aus C1D1 Pool rescuen
  semi_join(df.c1d1pat,Dynamicmut, by="Patmut")->DynamicCresc
  bind_rows(DynamicEresc,DynamicCresc)->Dynamicy
  
  #Matches und rescued zusammenführen
  bind_rows(Dynamicx,Dynamicy)->Dynamic
  
  #Bind Dynamic with other Pat Dynamics files
  bind_rows(Dynamicallpat,Dynamic)%>%
    unique()%>%
    filter(cf != "1")->Dynamicallpat
}


rm(C1D1EOTdynamic)+
  rm(df.pat)+
  rm(DynamicC)+
  rm(DynamicCresc)+
  rm(DynamicE)+
  rm(DynamicEresc)+
  rm(Dynamicmut)+
  rm(Dynamicx)+
  rm(Dynamicy)+
  rm(Dynamic)



#write to excel file
filename="C1D1EOT/Dynamicallpatients.xlsx"
#write.xlsx(Dynamicallpat,filename, sheetName = "Dynamic",append=TRUE)

#Create Dynamic plots for each EOT ID
#setwd to where images should be saved to!Muss noch ggeändert werdem, 
#sodass die plots in figures landen...
setwd("H:/Meine Ablage/C1D1EOT/Dynamicplot")
for (current_patient in EOT_ids)
{print(current_patient)
  Dynamicallpat%>%
    filter(Patient.ID == current_patient)->Dynamic
  ###Plot generieren
  Dynamic%>%
    ggplot(aes(y=TVAF, x=Visite, colour=Gene, group=Patmut))+
    geom_point(size=5,alpha=0.3)+
    geom_line(size=1)+
    theme_minimal()+
    scale_y_continuous(limits=c(0,0.125))+
    labs(title="Clone Dynamics")->p.C1D1EOTdynamicpat1
  ### create image file
  png(paste0(current_patient,"_C1D1EOT.png"),
      width=10,
      height=6,
      units="in",
      res=500,
      type="cairo")
  ### plot image to file
  print(p.C1D1EOTdynamicpat1)
  ### close file again
  dev.off()
}

rm(list=ls())
