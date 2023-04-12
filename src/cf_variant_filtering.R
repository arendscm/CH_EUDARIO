# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: from variantcalls to excel list of filtered results
#
# Input: seqdata
#
# Output: Excel list of filtered results
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

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')

##load hotspot data
load('data/interim/hotspots.RData')

##hotspot OvCa
#here we need data on mutational hotspots in ovca, e.g. from TCGA publication and all the others... I don't know if there are any?

######## Get Patient ids
source("src/ids.R")

########   FILTERING cf calls------------------------------------------------------------
df -> df.backup
df %>% filter(Material == "cf",
              Visite == "C1D1") -> df

##HRD genes
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3", "BRCA1", "BRCA2")
ovarian_cancer_genes <-c("TP53","NF1", "BRCA1", "BRCA2", "RB1","CDK12")
PARPi_actionable_genes <-c("ATM", "BRCA1", "BRCA2","BRIP1", "CDK12", "CHEK2", "PALB2")
ch_genes <- c("DNMT3A", "TET2" ,  "JAK2" ,  "ASXL1" , "SF3B1" , "SRSF2" , "TP53"  , "U2AF1" , "PPM1D" , "CBL"  ,  "IDH1"  , "IDH2"  , "BCOR"  , "BCORL1", "EZH2" ,  "RAD21" , "STAG2" , "CHEK2" , "GNAS"  , "GNB1"  , "ATM"   , "KRAS" ,  "NRAS",   "WT1" ,   "MYD88" ,
              "STAT3" , "BRCC3" , "CALR"  , "CEBPA" , "CSF3R" , "ETV6"  , "FLT3" ,  "GATA2" , "GATA1" , "KIT" ,   "MPL" ,   "NPM1" ,  "PTPN11" ,"RUNX1" , "SETBP1" ,"NF1"  ,  "PHF6")


# filter criteria
#functional criteria
df %>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  dplyr::select(mutID) %>% mutate(mutID = as.character(mutID)) -> mutID.func

## read count criteria
df %>%
  filter(readDepth > 100, TR2 > 19, TVAF > 0.005) %>%
  dplyr::select(mutID) -> mutID.count

## quality criteria
df %>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  dplyr::select(mutID)-> mutID.qual

## mutation call frequency criteria
df %>% 
  filter(AF<0.05) %>%    #filtert alle h?ufig in Datenbanken (=seq errors) gelisteten mutation raus
  filter((mutFreq < max(0.05*n.lane,5))&((p.binom<= -10)&med.vaf < 0.44))%>% #filtert nach Häufigkeit und binomialer Wahrscheinlichkeit
  dplyr::select(mutID) -> mutID.freq

# rescue ASXL1 dupG mutations <- this step is no longer needed, when we use p.binom 
#df %>%
#  filter(str_detect(Start,'32434638'))%>%filter(str_detect(AAChange,'ASXL1')) %>%
#  mutate(dev.med = ((TVAF - median(TVAF))/sd(TVAF))) %>% #calculate deviation from median in terms of standarddeviations
#  filter(dev.med > 1) %>%    #to be discussed
#  dplyr::select(mutID) -> mutID.asxl1

## rescue hotspots/mutations reported in cosmic data base for OVCa
df %>%
  filter(str_detect(cosmic92_coding,"ovary")) %>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 15) %>%
  filter(TVAF >0.001) %>%
  filter(!snp)%>%
  dplyr::select(mutID)-> mutID.cosmic


#filtering
##somatic variants
inner_join(mutID.func,mutID.count) %>% 
  inner_join(.,mutID.freq) %>% 
  inner_join(.,mutID.qual) %>% 
  inner_join(.,df) %>% 
  full_join(.,inner_join(df,mutID.cosmic))%>%
  filter(snp == FALSE) %>%
  #full_join(.,inner_join(df,mutID.tag.true))%>%
  filter(ExonicFunc != "synonymous SNV") %>% 
<<<<<<< HEAD
  filter(Gene %in% hrd_genes | Gene %in% ovarian_cancer_genes)%>%
=======
>>>>>>> 8c4a1c06d71fd61ed07d079c19fe7de6cbc61b1e
  group_by(Sample) %>% 
  mutate(n.mut.patient = n()) %>% 
  data.frame %>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding,"ovary"))%>%
  mutate(HRD = is.element(Gene,hrd_genes))%>%
  mutate(PARPi_actionable = is.element(Gene,PARPi_actionable_genes))%>%
  mutate(OvarianCancerGene = is.element(Gene,ovarian_cancer_genes))%>%
<<<<<<< HEAD
  mutate(CH_gene = is.element(Gene,ch_genes))-> df.filtered_cf 
=======
  mutate(CH_gene = is.element(Gene,ch_genes))%>%
  filter(!is.element(Sample.ID, failedSamples))-> df.filtered_cf
>>>>>>> 8c4a1c06d71fd61ed07d079c19fe7de6cbc61b1e

##hier kann man noch weiter filtern, hohe mutfreqs rausschmeißen, nur HRD Gene anschauen etc. und dann sollte es erstmal eine überschaubare menge an mutationen sein. Finetuning müssen wir dann noch schauen
##ich hab jetzt nur hrd und ovarian cancer genes, da wir über diese Tabelle ja die Patienten als hrd pos identififizieren... dazu gehören dann auch die pathogenen BRCA 1 & 2 germline mutations

<<<<<<< HEAD
rm(mutID.cosmic)+
rm(mutID.count)+
rm(mutID.freq)+
rm(mutID.func)+
rm(mutID.qual)+
=======
rm(mutID.count)+
rm(mutID.freq)+
rm(mutID.func)+
rm(mutID.cosmic)+
rm(mutID.tp53)+
rm(mutID.qual)+
rm(ids)+
>>>>>>> 8c4a1c06d71fd61ed07d079c19fe7de6cbc61b1e
rm(mm_hotspots)+
rm(ch_genes)+
rm(hrd_genes)+
rm(PARPi_actionable_genes)+
rm(ovarian_cancer_genes)+
rm(df.backup)+
rm(ids)+
rm(df)
rm(df.backup)

save.image("data/interim/seqdata_filtered_cf.RData")

<<<<<<< HEAD
filename <- paste("output/filtered_results_c1d1_cf_",Sys.Date(),".xlsx",sep="")
write.xlsx(df.filtered_cf,filename,sheetName = "filtered_results",append=TRUE)
=======
filename <- paste("output/filtered_results_c1d1_cf",Sys.Date(),".xlsx",sep="")
write.xlsx(df.filtered_cf_PARpi,filename,sheetName = "filtered_results",append=TRUE)



>>>>>>> 8c4a1c06d71fd61ed07d079c19fe7de6cbc61b1e
