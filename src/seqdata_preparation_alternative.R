# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: Variant call preprocessing
#
# Input: variantcalls as csv file, ids from ids.R
#
# Output: prepreprocessed variant list saved as either .csv or .RData
#
# ==============================================================================
########   Dependencies   #####
library(base)
library(dplyr)
library(xlsx)
library(reshape)
library(tidyr)
library(readxl)
library(reshape2)
library(stringr)
library(data.table)

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########   Data preparation:load variant calling table, Patient IDs, BRCA and tags-----------------------------------------------------------------
data1 <- read.table('data/raw/Run1/variantcalls_P1346.csv',
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)


data2 <- read.table('data/raw/Run2/variantcalls_P1519.csv',
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)


data3 <- read.table('data/raw/Run3/variantcalls_P1803.csv',
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#load run and lane data
lanes <- read_excel("data/external/Sample-Run-Assignment.xlsx")%>%
  mutate(lane = paste(Run,Lane,sep="_"))


# load tags 
#tags <- read.table('data/external/tags_run1_run2.csv',
#                   header = TRUE, sep = ";", stringsAsFactors = FALSE)
load("data/interim/tags.RDATA")
tags%>% 
  filter(!is.element(mutID,c("1-G7_chr17_7674872_7674872_T_C","3-C1_chr17_60663388_60663388_C_T","3-F7_chr17_7674220_7674220_C_T")))-> tags

##Patient ID table that identifies Sample IDs with Patient ID and timepoints
source("src/material_table.R")


########   Data Preparation: Restructure data for filtering---------------------------------------------------------------------
##relevant variables
variables <- c("Sample", "Patient.ID", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "GeneDetail", "ExonicFunc", "AAChange", "cytoBand","readDepth", "TR1", "TR1_plus", "TR1_minus", "TR2", "TR2_plus", "TR2_minus", "TVAF", "AF", "avsnp150", "cosmic92_coding")
discard_variables <- c("CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG", "SIFT_score", "SIFT_converted_rankscore", "SIFT_pred", "LRT_score", "LRT_converted_rankscore", "LRT_pred", "MutationTaster_score", "MutationTaster_converted_rankscore", "MutationTaster_pred", "MutationAssessor_score", "MutationAssessor_score_rankscore", "MutationAssessor_pred", "FATHMM_score", "FATHMM_converted_rankscore", "FATHMM_pred", "PROVEAN_score", "PROVEAN_converted_rankscore", "PROVEAN_pred", "MetaSVM_score", "MetaSVM_rankscore", "MetaSVM_pred", "MetaLR_score", "MetaLR_rankscore", "MetaLR_pred", "M.CAP_score", "M.CAP_rankscore", "M.CAP_pred", "MutPred_score", "MutPred_rankscore", "fathmm.MKL_coding_score", "fathmm.MKL_coding_rankscore", "fathmm.MKL_coding_pred", "Eigen_coding_or_noncoding", "Eigen.raw", "Eigen.PC.raw", "GenoCanyon_score", "GenoCanyon_score_rankscore", "integrated_fitCons_score", "integrated_fitCons_score_rankscore", "integrated_confidence_value", "GERPpp_RS", "GERP.._RS_rankscore", "phyloP100way_vertebrate", "phyloP100way_vertebrate_rankscore", "phyloP20way_mammalian", "phyloP20way_mammalian_rankscore", "phastCons100way_vertebrate", "phastCons100way_vertebrate_rankscore", "phastCons20way_mammalian", "phastCons20way_mammalian_rankscore", "SiPhy_29way_logOdds", "SiPhy_29way_logOdds_rankscore", "Interpro_domain", "GTEx_V6p_gene", "GTEx_V6p_tissue", "fwd.Primer", "rev_Primer", "Status", "Note", "AmpliconRange", "AmpliconSize", "InsertRange", "InsertSize", "InsertSeq", "offsetL", "offsetR")

#hotspots
hotspots <- c("chr2_25234374_25234374",     #DNMT3A R882C
              "chr2_25234373_25234373",     #DNMT3A R882H
              "chr2_208248389_208248389",   #IDH1 R132
              "chr1_1815790_1815790",       #GNB1 K57
              "chr1_1815789_1815789",       #GNB1 K57
              "chr12_25245350_25245350",    #KRAS G12D
              "chr2_197402636_197402636",   #SF3B1 K666
              "chr2 197402635 197402635",   #SF3B1 K666
              "chr2_197402110_197402110",   #SF3B1 K700
              "chr17_76736877_76736877",    #SRSF2 P95
              "chr9_5073770_5073770",       #JAK2 V617
              "chr20_58909366_58909366"     #GNAS R844
)

#fuse inputdata and join with run and lane info
data <- full_join(data1%>%dplyr::select(-discard_variables),data2%>%dplyr::select(-discard_variables))%>%
  full_join(.,data3%>%dplyr::select(-discard_variables))%>%
  full_join(.,lanes)

##calculate mutation frequencies lane-wise
data %>% 
  mutate(Sample_orig = Sample)%>%
  left_join(.,ids,by = "Sample_orig") %>% ##fuse with SampleID Table
  mutate(AF = replace(AF, AF == ".", 0)) %>% 
  mutate(AF = as.numeric(AF)) %>%
  mutate(position = paste(Chr,Start,End,Ref,Alt,sep="_"))%>%  #mutation position with basechange
  mutate(position2 = paste(Chr,Start,End,sep="_"))%>% #mutation position without basechange
  mutate(mutID = paste(Sample_orig,position,sep="_")) %>%    #mutation id (here we use Sample_orig, because it is unique)
  mutate(freqID = paste(position,AAChange,sep=":")) %>%   #id to calculate mutation frequency
  mutate(snp = ((TVAF > 0.4 & TVAF < 0.6)|(TVAF > 0.9))&(AF>0.001))%>% #SNP classification
  group_by(lane)%>%
  mutate(mutFreq = as.vector(table(freqID)[freqID])) %>%  #mutation frequency
    data.frame -> tmp

#clear memory space
rm(list(data1,data2,data3,data))

#modify/add columns needed for filtering
df <- data.frame(tmp) %>% 
  full_join(.,tags) %>% #join with list of tags from manual inspection in igv
  mutate(tag = as.factor(tag)) %>%
  unique %>%
  group_by(Patient.ID,position) %>% 
  mutate(n.mut = n()) %>% #add column stating whether mutation is present in several samples from the same Patient.ID
  mutate(Genomic_Coordinate_hg38 = paste(Chr,paste("g",Start,sep = "."),paste(Ref,Alt,sep =">"),sep=":"))%>% #Genomic coordinate to match with BRCA_exchange database
  data.frame %>%
  group_by(Patient.ID,position,Material) %>%
  mutate(n.material = n())%>% #to see number of times this mutation occurs within the same material (= at different timepoints)
  data.frame

  as.data.frame(df %>% 
                group_by(freqID) %>% 
                dplyr::summarise(mutID,n(),TVAF,median(TVAF),mean(TVAF),sd(TVAF))) %>% #for every mutation position (as determined by freqID), calculate median VAF, mean VAF and standard deviation. This will be relevant for positions with high mutation frequency
  dplyr::select(mutID,TVAF,'median(TVAF)','sd(TVAF)','mean(TVAF)') %>% 
  dplyr::rename(med.vaf = 'median(TVAF)',sd.vaf = 'sd(TVAF)',mean.vaf = 'mean(TVAF)') %>%
  left_join(df,.) %>% 
  mutate(dev.med = (TVAF-med.vaf)/sd.vaf) %>% #tells us by how many standard deviations the VAF at the current position deviates from the median
  mutate(p.binom = ifelse(mutFreq > 9,log10(1-pbinom(TR2,readDepth,med.vaf)),-Inf))%>%   #p.binom is a measure for the likelihood to get these number of variant reads TR2 (or more), assuming a Bernoulli experiment with p=med.vaf
  filter(!is.na(TR2),!is.na(readDepth))%>% 
  group_by(Patient.ID,position) %>% 
  mutate(serial.mut = n()) %>%
  data.frame %>%
  mutate(hotspot = is.element(position2,hotspots)) %>% #hotspots (defined above)
  mutate(igv = paste(Chr,Start,sep=":"))%>%
  mutate(Patmut=paste(Patient.ID,position,sep="_"))%>%
    left_join(.,df.material,by="Patient.ID")->df

#save as csv in interim data folder
write.csv(df%>%filter(!is.na(Patient.ID)),'data/interim/mutationcalls.csv')
write.csv(df%>%filter(is.na(Patient.ID)),'data/interim/newsample_calls.csv')

##Save RData for further use
  tempdata <-ls()
  rm(list=tempdata[tempdata != "df"])
  rm(tempdata)
  
save.image('data/interim/seqdata.RData')
  