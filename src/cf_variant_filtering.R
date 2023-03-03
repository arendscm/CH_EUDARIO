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
  filter(readDepth > 100, TVAF > 0.01, TR2>19) %>%
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

## rescue hotspots/mutations reported in cosmic data base for OVCa
df %>%
  filter(str_detect(cosmic92_coding,"ovary")) %>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 9) %>%
  filter(TVAF >0.001) %>%
  filter(p.binom< -10)%>%
  filter(mutFreq < 0.1*n.lane)%>%
  dplyr::select(mutID)-> mutID.cosmic

<<<<<<< HEAD
=======
## rescue TP53 mutations OVCa
df %>%
  filter(Gene == "TP53") %>%
  filter(ExonicFunc != "synonymous SNV")%>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 14) %>%
  filter(TVAF >0.001) %>%
  filter(p.binom< -10)%>%
  filter(mutFreq < 0.1*n.lane)%>%
  dplyr::select(mutID)-> mutID.tp53
>>>>>>> 49d04a274793327386bc0cf46e16aeedcba7869a

#filtering
##somatic variants
inner_join(mutID.func,mutID.count) %>% 
  inner_join(.,mutID.freq) %>% 
  inner_join(.,mutID.qual) %>% 
  inner_join(.,df) %>% 
  full_join(.,inner_join(df,mutID.cosmic))%>%
  full_join(.,inner_join(df,mutID.tp53))%>%
  filter(snp == FALSE&!(TVAF > 0.4 & n.visite>1)) %>% #filters out SNPs and germline mutations
  #full_join(.,inner_join(df,mutID.tag.true))%>%
  filter(ExonicFunc != "synonymous SNV") %>% 
  filter(Gene %in% hrd_genes | Gene %in% ovarian_cancer_genes)%>%
  group_by(Sample) %>% 
  mutate(n.mut.patient = n()) %>% 
  data.frame %>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding,"ovary"))%>%
  mutate(HRD = is.element(Gene,hrd_genes))%>%
  mutate(PARPi_actionable = is.element(Gene,PARPi_actionable_genes))%>%
  mutate(OvarianCancerGene = is.element(Gene,ovarian_cancer_genes))%>%
  mutate(CH_gene = is.element(Gene,ch_genes))-> df.filtered_cf 

##hier kann man noch weiter filtern, hohe mutfreqs rausschmeißen, nur HRD Gene anschauen etc. und dann sollte es erstmal eine überschaubare menge an mutationen sein. Finetuning müssen wir dann noch schauen
##ich hab jetzt nur hrd und ovarian cancer genes, da wir über diese Tabelle ja die Patienten als hrd pos identififizieren... dazu gehören dann auch die pathogenen BRCA 1 & 2 germline mutations



rm(mutID.CHIP)+
rm(mutID.CHIP.qual)+
rm(mutID.count)+
rm(mutID.freq)+
rm(mutID.func)+
rm(mutID.hotspots)+
rm(mutID.qual)+
rm(mutID.tag.true)
rm(ids)+
rm(mm_hotspots)+
rm(tempdata)+
rm(df)

save.image("data/interim/seqdata_filtered_cf.RData")

filename <- paste("output/filtered_results_c1d1_cf_",Sys.Date(),".xlsx",sep="")
write.xlsx(df.filtered_cf,filename,sheetName = "filtered_results",append=TRUE)

#####BRCA somatic mutations
df %>%
  filter(Gene == "BRCA1"|Gene == "BRCA2") %>%
  filter(ExonicFunc != "synonymous SNV")%>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 19) %>%
  filter(TVAF >0.01) %>%
  filter(p.binom< -10)%>%
  filter(mutFreq < 0.1*n.lane)%>%
  filter(!snp)


