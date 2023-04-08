# ______________________________________________________________________________
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: determine mutation signatures in cfDNA mutations and wb mutations. 
#
# Input: df from data/interim/seqdata.RData, ...
#
# Output: plots...
#
# ______________________________________________________________________________
#####  Dependencies   #####
library(base)
library(dplyr)
library(stringr)
library(reshape)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(viridis)
library(ggpubr)
#####  Data preparation ####
########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')
load('data/interim/seqdata_filtered.RData')

######## Get Patient ids
source("src/ids.R")
source("src/material_table.R")

######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")

##interesting gene groups
variables <- c("Patient.ID","Sample_orig","mutID","position","Sample", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "GeneDetail", "ExonicFunc", "AAChange", "cytoBand","readDepth", "TR1", "TR1_plus", "TR1_minus", "TR2", "TR2_plus", "TR2_minus", "TVAF", "AF", "avsnp150","cosmic92_coding","snp","mutFreq","p.binom","n.mut","n.material","sum_cf","sum_wb","Material","tag", "Patmut")
ch_genes_without_HRD <- c("DNMT3A","TET2","ASXL1","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5")
tp53_genes <- c("TP53")
ppm1d_genes <- c("PPM1D")
brca_genes <- c("BRCA1","BRCA2")
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")
failedSamples <-c('OvCA_44_C1D1_cf','OvCA_45_C1D1_cf','OvCA_46_C1D1_cf','OvCA_48_C1D1_cf','OvCA_50_C1D1_cf','OvCA_54_C1D1_cf','OvCA_93_C1D1_cf',
                  'OvCA_11_C1D1_cf','OvCA_40_C1D1_cf','OvCA_53_C1D1_cf','OvCA_65_C1D1_cf')
Categories<-c('CH','HRD','other','TP53')

#dataframe with mutation calls from cfDNA             
df %>% 
  filter(Material=="cf") %>% 
  filter(Visite == "C1D1")%>%
  filter(is.na(replicate))%>%
  filter(!is.element(Sample.ID, failedSamples))%>%
  dplyr::select(all_of(variables)) %>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))-> df.cf

#data frame with mutation calls from WB samples that have matched cfDNA samples
df %>% 
  filter(is.element(Patient.ID,df.cf$Patient.ID)) %>% 
  filter(is.na(replicate))%>%
  filter(Material=="wb") %>% 
  filter(Visite == "C1D1")%>%
  filter(!is.element(Sample.ID, failedSamples))%>%
  dplyr::select(all_of(variables)) %>%
  mutate(cfID=paste(Patient.ID,position,sep="_")) -> df.cf_wb

#####  Determine trinucleotide context of mutations cf vs wb


BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library(deconstructSigs)
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(!is.element(Sample.x,failedSamples))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes_without_HRD),"CH",
                       ifelse(is.element(Gene.x,tp53_genes),"TP53",
                              ifelse(is.element(Gene.x,hrd_genes),"HRD",
                                     ifelse(is.element(Gene.x,brca_genes),"BRCA",
                                            ifelse(is.element(Gene.x,ppm1d_genes),"PPM1D","other"))))))%>%
  #filter(gene != "other") %>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding.x,"ovary")) %>%
  filter(p.binom.x <= -Inf) %>%
  filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  filter(ExonicFunc.x != "synonymous SNV")%>%
  filter(AF.x<0.1)%>%
  filter(snp.x==FALSE)%>%
  filter(TVAF.x>0.01|TVAF.y>0.01)%>%
  filter(TR2.y > 19|TR2.x>19)%>%
  filter(TVAF.x<0.35&TVAF.y<0.35)%>% #no germline variants
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  mutate(chr = Chr.x,
         sample.id = compartment,
         pos = Start.x,
         ref = Ref.x,
         alt= Alt.x)%>%
  filter(ExonicFunc.x=="nonsynonymous SNV")%>%
  dplyr::select(c("sample.id","chr","pos","ref","alt"))%>%
  mut.to.sigs.input(mut.ref = ., 
                    sample.id = "sample.id", 
                    chr = "chr", 
                    pos = "pos", 
                    ref = "ref", 
                    alt = "alt")

