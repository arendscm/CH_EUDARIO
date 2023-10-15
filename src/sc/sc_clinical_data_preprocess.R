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
library(fastDummies)


########   Load IDs    ########

source("src/genegroup_definitions.R")

########## Load mutation data with clinical data ##########
load("data/interim/newsamples_filtered.RData")

SC_registry <- read.xlsx("data/external/SC_registry.xlsx",sheetIndex=1,header=TRUE)


########   Load clinical data    ########
df.clin_orig <- data.frame(read.xlsx("data/external/SC_clindata.xlsx",1,header=TRUE))

##load BRCA germline status
load("data/interim/brca_sc.Rdata")

#irrelevant variables
irrel_var <- c("therapy.prior.d1","Bemerkung","Cancer.in.history","HRD")

##preprocess df.clin
df.clin_sc <- df.clin_orig %>% 
  dplyr::select(-irrel_var)%>%
  mutate(no_prev_lines_binom = factor(ifelse(priorLines == 1,"1",">1"),levels = c("1",">1")),
         Duration_PARPi_level = factor(ifelse(Duration_PriorPARPi == 0,"0",
                                       ifelse(Duration_PriorPARPi<10,"<10 months",">10 months")),levels=c("0","<10 months",">10 months")))%>%
  left_join(.,id.brca_germline_SC, by="Patient.ID") %>% #add BRCA germline status
  left_join(.,id.hrd_germline_SC, by="Patient.ID") %>% #add HRD germline status
  mutate(HRD_germline = sign(brca1_germline+brca2_germline+hrd_germline),
         brca_germline = sign(brca1_germline+brca2_germline))

########## Join mutation data with clinical data ##########
df.mut <- df.filtered.newsamples%>% 
  filter(TVAF >= 0.01)%>%
  filter(tag=="true")%>%
  dplyr::select("Patient.ID","Sample", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "ExonicFunc", "AAChange","TR1","TR2","TVAF","cosmic92_coding") %>%
  mutate(CH = is.element(Gene,typical_ch_genes),
         HRD = is.element(Gene,hrd_genes),
         DDR = is.element(Gene,ddr_genes))%>%
  group_by(Patient.ID,CH)%>%
  mutate(maxVAF_CH = max(TVAF))%>%
  data.frame%>%
  mutate(maxVAF_CH = ifelse(CH ==1, maxVAF_CH,0))%>%
  group_by(Patient.ID)%>%
  mutate(maxVAF_CH = max(maxVAF_CH))%>%
  data.frame

#create dummies
df.mut %>% 
  mutate(Gene = ifelse(Gene == "U2AF1;U2AF1L5","U2AF1",Gene))%>%
  dummy_cols(., select_columns = "Gene")%>%
  group_by(Patient.ID)%>%
  summarise(nom = n(), 
            nom_CH = sum(CH),
            DNMT3A = sum(Gene_DNMT3A),
            TET2 = sum(Gene_TET2),
            ASXL1 = sum(Gene_ASXL1),
            TP53 = sum(Gene_TP53),
            PPM1D = sum(Gene_PPM1D),
            CHEK2 = sum(Gene_CHEK2),
            maxVAF = max(TVAF),
            maxVAF_CH = max(maxVAF_CH),
            CH = sign(sum(CH)),
            HRD = sign(sum(HRD)),
            DDR = sign(sum(DDR)))%>%
  data.frame%>%
  mutate(DTA = sign(DNMT3A + TET2 + ASXL1),
         tp53ppm1d = sign(TP53+PPM1D),
         any_clone = 1,
         CHIP = ifelse(maxVAF_CH >= 0.02,1,0),
         CH5 = ifelse(maxVAF_CH >= 0.05,1,0),
         CH10 = ifelse(maxVAF_CH >= 0.1,1,0),
         mult.mut = ifelse(nom_CH > 1,2,ifelse(nom == 0,0,1)),
         nond3a.mut = sign(nom_CH - DNMT3A)
  ) -> df.mut.dummy

##Fuse mutation data with clinical data
df <- df.clin_sc %>% 
  full_join(.,df.mut.dummy,by="Patient.ID") %>% 
  mutate_at(.vars = vars(c("DNMT3A", "TET2", "ASXL1","PPM1D","TP53","CHEK2","tp53ppm1d","DTA","DDR","CH","any_clone","HRD","CHIP","CH5","CH10","nond3a.mut")),
            .funs = list(~as.factor(sign(ifelse(is.na(.),0,.))))) %>%  
  mutate_at(.vars = vars(c("nom","nom_CH","maxVAF","maxVAF_CH","mult.mut")),
            .funs = list(~ifelse(is.na(.),0,.))) %>%
  mutate(CH_category = cut(maxVAF_CH,breaks=c(-Inf,0.01,0.1,Inf),labels=c("<1%","1-10%",">10%")))

##Save RData for further use
tempdata <-ls()
rm(list=tempdata[tempdata != "df"])
rm(tempdata)

df -> df.clin_sc
rm(df)

save.image('data/interim/sc_clin.RData')

