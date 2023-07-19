########   MAFTOOLS - Oncoplot with clinical annotation#### 

library(base)
library(dplyr)
library(stringr)
library(reshape)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(reshape)
library(ggpubr)
library(maftools)
library(RColorBrewer)
library(ggsci)

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata_filtered.RData')
load('data/interim/clin.RData')

######## Get Patient ids
source("src/material_table.R")

######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")
source("src/genegroup_definitions.R")


###turn df.filtered into MAF compatible format

maf<-full_join(ids %>% 
                  filter(firstTimepoint_wb==1,is.na(replicate)) %>% 
                  mutate(Sample_ID=Patient.ID, Tumor_Sample_Barcode=Patient.ID) %>% 
                  dplyr::select(Sample_ID,Tumor_Sample_Barcode),makeMAF(df.filtered.c1d1%>% 
                                                                          mutate(Gene=ifelse(Gene=="U2AF1;U2AF1L5","U2AF1",Gene))%>%
                                                                          filter(tag=="true",TVAF >= 0.01))) 
clin<- df.clin %>% 
              mutate(Sample_ID = Patient.ID,
                     Tumor_Sample_Barcode = Patient.ID,
                     Age = Age_TreatmentStartEUDARIO,
                     Arm = ifelse(Arm==" A","A",
                                  ifelse(Arm==" B","B","C")))%>%
              dplyr::select(Sample_ID,Tumor_Sample_Barcode,Age,Arm,PriorPARPi,Number_PreviousLines )

df.maf <- read.maf(maf,clinicalData = clin)

vc_cols = pal_npg("nrc")(8)
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

Arm_colors=pal_npg("nrc")(3)
names(Arm_colors)=c('A','B','C')

PARPi_colors = pal_npg("nrc")(3)
names(PARPi_colors)=c( "Yes", "No", NA )

Lines_colors=scales::seq_gradient_pal(high = "#E64B35FF", low = "#4DBBD5FF",space="Lab")(seq(0,1,length.out=5))
names(Lines_colors)=c(1,2,3,4,5)


#This is the final colored list. Names of the list elements should match those in clinicalFeatures arguments 
anno_cols = list(Arm = Arm_colors, PriorPARPi = PARPi_colors, Number_PreviousLines = Lines_colors)

pw1<-data.frame(ch_genes_without_HRD,"CH")
names(pw1) <- c("Genes","Pathway")
pw2<-data.frame(hrd_genes,"HRD")
names(pw2) <- c("Genes","Pathway")
Group <- rbind(pw1,pw2)


png("output/figures/oncoplot.png",width=12, height=8,units="in",res=500,type="cairo")

oncoplot(df.maf,
         top = 22,
         colors = vc_cols,
         drawColBar = TRUE,
         topBarData = df.clin %>% mutate(Tumor_Sample_Barcode = Patient.ID, 'No. of mutations' = nom) %>% dplyr::select(Tumor_Sample_Barcode,'No. of mutations'),
         removeNonMutated = FALSE,
         clinicalFeatures = c('PriorPARPi','Number_PreviousLines','Arm'),
         annotationColor = anno_cols,
         pathways = Group,
         sortByAnnotation = TRUE,
         annotationOrder = c("Yes","No","NA"),
         anno_height = 1.5) 

dev.off()

##somatic interactions
somaticInteractions(df.maf, top=5)