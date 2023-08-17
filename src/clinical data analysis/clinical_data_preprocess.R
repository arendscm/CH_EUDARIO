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
source("src/material_table.R")
source("src/genegroup_definitions.R")

########## Load mutation data with clinical data ##########
load("data/interim/seqdata_filtered.RData")

########   Load clinical data    ########
df.clin_orig <- data.frame(read.xlsx("data/external/clindata_modified.xlsx",1,header=TRUE))
df.ae <- data.frame(read.table("data/external/ae_modified.csv",header=TRUE,sep=";"))
df.bc <- data.frame(read.table("data/external/bc_modified.csv",header=TRUE,sep=";"))
df.response <- data.frame(read.table("data/external/response_modified.csv",header=TRUE,sep=";"))

##load BRCA germline status
load("data/interim/brca.Rdata")
load("data/interim/brca_somatic.Rdata")

#irrelevant variables
irrel_var <- c("Date_subsequent.progression","AdditionalCancerTherapy","ORR")

##preprocess df.clin
df.clin <- df.clin_orig %>% 
  dplyr::select(-irrel_var)%>%
  mutate(ECOG = factor(ECOG_SCR),
         ECOG_binom = factor(ifelse(ECOG == 0,0,1)),
         TumorBurden_baseline = as.numeric(TumorBurden_baseline),
         LVEF_C1D1 = as.numeric(LVEF_C1D1),
         No_Platinum_lines_binom = factor(ifelse(Number_PreviousPlatinumLines>1,">1","1"),levels=c("1",">1")),
         no_prev_lines_binom = factor(ifelse(Number_PreviousLines == 1,"1",">1"),levels = c("1",">1")),
         Duration_PriorPARPi = ifelse(is.na(Duration_PriorPARPi),0,Duration_PriorPARPi),
         Duration_PARPi_level = factor(ifelse(Duration_PriorPARPi == 0,"0",
                                       ifelse(Duration_PriorPARPi<10,"<10 months",">10 months")),levels=c("0","<10 months",">10 months")),
         response_binom = ifelse(Response_best == "NAP",NA,
                                 ifelse(Response_best =="CR"|Response_best =="PR","CR or PR","SD or PD")),
         response_binom2 = ifelse(Response_best == "NAP",NA,
                                 ifelse(Response_best =="CR"|Response_best =="PR"|Response_best =="SD","CR, PR or SD","PD")),
         response_binom3 = ifelse(Response_best == "NAP",NA,
                           ifelse(Response_best =="CR","CR","no CR")))%>%
  mutate(ae_haematotox = is.element(Patient.ID,(df.ae %>%
                                                  filter(is.element(AE_consensus_term,c("neutropenia","thrombocytopenia","anemia")))%>%
                                                  dplyr::select(patient_code)%>% 
                                                  unique)$patient_code),
         ae_neutropenia = is.element(Patient.ID,(df.ae %>%
                                                   filter(is.element(AE_consensus_term,c("neutropenia")))%>%
                                                   dplyr::select(patient_code)%>% 
                                                   unique)$patient_code),
         ae_thrombocytopenia = is.element(Patient.ID,(df.ae %>%
                                                   filter(is.element(AE_consensus_term,c("thrombocytopenia")))%>%
                                                   dplyr::select(patient_code)%>% 
                                                   unique)$patient_code),
         ae_anemia = is.element(Patient.ID,(df.ae %>%
                                                        filter(is.element(AE_consensus_term,c("anemia")))%>%
                                                        dplyr::select(patient_code)%>% 
                                                        unique)$patient_code),
         ae_haematotox_severe = is.element(Patient.ID,(df.ae %>%
                                                         filter(is.element(AE_consensus_term,c("neutropenia","thrombocytopenia","anemia"))&AE_18_severity>2)%>%
                                                         dplyr::select(patient_code)%>% 
                                                         unique)$patient_code),
         ae_infection = is.element(Patient.ID,(df.ae %>%
                                                         filter(is.element(AE_consensus_term,c("infection")))%>%
                                                         dplyr::select(patient_code)%>% 
                                                         unique)$patient_code),
         ae_infection_severe = is.element(Patient.ID,(df.ae %>%
                                                 filter(is.element(AE_consensus_term,c("infection"))&AE_18_severity>2)%>%
                                                 dplyr::select(patient_code)%>% 
                                                 unique)$patient_code),
         ae_allergic_reaction = is.element(Patient.ID,(df.ae %>%
                                                         filter(is.element(AE_consensus_term,c("allergic reaction")))%>%
                                                         dplyr::select(patient_code)%>% 
                                                         unique)$patient_code),
         ae_allergic_reaction_severe = is.element(Patient.ID,(df.ae %>%
                                                         filter(is.element(AE_consensus_term,c("allergic reaction"))&AE_18_severity>2)%>%
                                                         dplyr::select(patient_code)%>% 
                                                         unique)$patient_code),
         ae_kidney_failure = is.element(Patient.ID,(df.ae %>%
                                                                filter(is.element(AE_consensus_term,c("increased creatinine","renal insufficiency","kidney insufficiency")))%>%
                                                                dplyr::select(patient_code)%>% 
                                                                unique)$patient_code),
         ae_transaminases = is.element(Patient.ID,(df.ae %>%
                                                      filter(is.element(AE_consensus_term,c("increased transaminases")))%>%
                                                      dplyr::select(patient_code)%>% 
                                                      unique)$patient_code),
         ae_bleeding = is.element(Patient.ID,(df.ae %>%
                                                filter(is.element(AE_consensus_term,c("bleeding")))%>%
                                                dplyr::select(patient_code)%>% 
                                                unique)$patient_code)
         )%>%
  filter(is.element(Patient.ID,ids %>% filter(firstTimepoint_wb==1) %>%.$Patient.ID)) %>% #only patients with available seq data
  left_join(.,id.brca_germline, by="Patient.ID") %>% #add BRCA germline status
  left_join(.,id.brca_somatic, by="Patient.ID") %>% #add BRCA somatic muts
  left_join(.,id.hrd_germline, by="Patient.ID")%>%
  mutate(HRD_germline = sign(brca1_germline+brca2_germline+hrd_germline),
         brca_germline = sign(brca1_germline+brca2_germline),
         BRCA1 = ifelse(BRCA1==" ",-1,
                        ifelse(BRCA1==" wt",0,1)),
         BRCA2 = ifelse(BRCA2==" ",-1,
                        ifelse(BRCA2==" wt",0,1)))%>%
  mutate(BRCA_status = ifelse(brca_germline ==1, "mut",
                              ifelse(BRCA1==1|BRCA2==1,"mut",
                                     ifelse(BRCA1==0&BRCA2==0,"wt","unknown"))),
         BRCA_binom = ifelse(BRCA_status=="mut","mut","wt or unknown"))%>%
  left_join(.,df.bc%>%
              mutate(Patient.ID = as.character(patient_code))%>%
              filter(cycle_day==" C1D1")%>%
              dplyr::select(Patient.ID,hemoglobin,thr,wbc,CA125))

########## Join mutation data with clinical data ##########
df.mut <- df.filtered.c1d1%>% 
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
            SF3B1 = sum(Gene_SF3B1),
            U2AF1 = sum(Gene_U2AF1),
            TP53 = sum(Gene_TP53),
            PPM1D = sum(Gene_PPM1D),
            CHEK2 = sum(Gene_CHEK2),
            ATM = sum(Gene_ATM),
            maxVAF = max(TVAF),
            maxVAF_CH = max(maxVAF_CH),
            CH = sign(sum(CH)),
            HRD = sign(sum(HRD)),
            DDR = sign(sum(DDR)))%>%
  data.frame%>%
  mutate(DTA = sign(DNMT3A + TET2 + ASXL1),
         tp53ppm1d = sign(TP53+PPM1D),
         splice_mut = sign(SF3B1+U2AF1),
         any_clone = 1,
         CHIP = ifelse(maxVAF_CH >= 0.02,1,0),
         CH5 = ifelse(maxVAF_CH >= 0.05,1,0),
         CH10 = ifelse(maxVAF_CH >= 0.1,1,0),
         mult.mut = ifelse(nom_CH > 1,2,ifelse(nom == 0,0,1)),
         nond3a.mut = sign(nom_CH - DNMT3A)
  ) -> df.mut.dummy

##Fuse mutation data with clinical data
df <- df.clin %>% 
  full_join(.,df.mut.dummy,by="Patient.ID") %>% 
  mutate_at(.vars = vars(c("DNMT3A", "TET2", "ASXL1","SF3B1","U2AF1","PPM1D","TP53","CHEK2","ATM","tp53ppm1d","DTA","DDR","splice_mut","CH","any_clone","HRD","CHIP","CH5","CH10","nond3a.mut")),
            .funs = list(~as.factor(sign(ifelse(is.na(.),0,.))))) %>%  
  mutate_at(.vars = vars(c("nom","nom_CH","maxVAF","maxVAF_CH","mult.mut")),
            .funs = list(~ifelse(is.na(.),0,.))) %>%
  mutate(CH_category = cut(maxVAF_CH,breaks=c(-Inf,0.01,0.1,Inf),labels=c("<1%","1-10%",">10%")))

##Save RData for further use
tempdata <-ls()
rm(list=tempdata[tempdata != "df"])
rm(tempdata)

df -> df.clin
rm(df)

save.image('data/interim/clin.RData')

