# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: mutational spectrum in patients with and without prior PARPi
#
# Input: preprocessed clinical data and mutational data
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
library(tableone)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggmosaic)

########   Load IDs    ########
source("src/material_table.R")
source("src/global_functions_themes.R")
source("src/genegroup_definitions.R")

########## Load clinical data  ##########
load("data/interim/clin.RData")
load("data/interim/seqdata_filtered.RData")


######## Plot with mut spectrum in pts with and without prior PARPi

df.clin$PriorPARPi %>% table %>% data.frame -> n.parp
names(n.parp) <- c("PriorPARPi","n.parp")

genes=c("DNMT3A","TET2","ASXL1","PPM1D","TP53","CHEK2","ATM")

df.filtered.c1d1 %>% 
  left_join(.,df.clin %>% 
              group_by(PriorPARPi) %>% 
              mutate(n.parp = n()) %>% 
              data.frame %>% 
              dplyr::select(Patient.ID, PriorPARPi,n.parp)%>%
              mutate(PriorPARPi = factor(PriorPARPi,levels=c("Yes","No"))))%>%
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  filter(Gene %in% genes)%>%
  dplyr::select(Sample, Gene, PriorPARPi,n.parp) %>% 
  data.frame %>% 
  unique %>% 
  dplyr::select(Gene,PriorPARPi) %>% 
  table %>% 
  data.frame %>% 
  left_join(.,n.parp)%>%
  mutate(prev = Freq/n.parp) %>% 
  arrange(prev) %>%
  ggplot(aes(x=reorder(Gene, Freq), y=prev, fill=PriorPARPi)) +
  geom_bar(stat="identity", width=0.6, position=position_dodge())+
  xlab("")+
  scale_y_continuous(labels = percent,limits=c(0,0.6), position = "right")+
  ylab("Gene Mutation Prevalence [%]") +
  my_theme()+
  scale_fill_npg(name="Prior PARPi therapy")+
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  coord_flip() -> p.mutprev_parpi
p.mutprev_parpi

png("output/figures/mutprev_priorparpi.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev_parpi
dev.off()

####### mut spect in pats with differing no of prior platinum lines

df.clin$No_Platinum_lines_binom %>% table %>% data.frame -> n.plat
names(n.plat) <- c("No_Platinum_lines_binom","n.plat")

genes=c("DNMT3A","TET2","ASXL1","PPM1D","TP53","CHEK2","ATM")

df.filtered.c1d1 %>% left_join(.,df.clin  %>% dplyr::select(Patient.ID, No_Platinum_lines_binom))%>%
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  filter(Gene %in% genes)%>%
  dplyr::select(Sample, Gene, No_Platinum_lines_binom) %>% 
  data.frame %>% 
  unique %>% 
  dplyr::select(Gene,No_Platinum_lines_binom) %>% 
  table %>% 
  data.frame %>% 
  left_join(.,n.plat)%>%
  mutate(prev = Freq/n.plat) %>% 
  arrange(prev) %>%
  ggplot(aes(x=reorder(Gene, Freq), y=prev, fill=No_Platinum_lines_binom)) +
  geom_bar(stat="identity", width=0.6, position=position_dodge())+
  xlab("")+
  scale_y_continuous(labels = percent,limits=c(0,0.6), position = "right")+
  ylab("Gene Mutation Prevalence [%]") +
  my_theme()+
  scale_fill_npg()+
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  coord_flip() -> p.mutprev_plat
p.mutprev_plat
  

png("output/figures/mutprev_priorplatinum.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev_plat
dev.off()

#####determinants of ch

my_vars_baseline=c("Number_PreviousLines",
                   "Number_PreviousPlatinumLines",
                   "no_prev_lines_binom",
                   "No_Platinum_lines_binom",
                   "Type_PreviousTherapy",
                   "PriorPARPi",
                   "Duration_PriorPARPi",
                   "Duration_PARPi_level")

cat_vars_baseline=c("Type_PreviousTherapy",
                    "PriorPARPi",
                    "no_prev_lines_binom",
                    "Duration_PARPi_level",
                    "No_Platinum_lines_binom")

cont_vars_baseline = setdiff(my_vars_baseline,cat_vars_baseline)

df.clin %>% CreateTableOne(strata = "CH",
                      vars=c(my_vars_baseline),
                      factorVars = cat_vars_baseline,
                      includeNA=FALSE,
                      #addOverall = TRUE,
                      data=.) %>% 
  print(., 
        #nonnormal=cont_vars_baseline,
        exact=cat_vars_baseline,
        missing=TRUE,
        showAllLevels=TRUE,
        quote=FALSE)


##### compare number of mutations as a function of prior PARPi

df.clin %>% 
  dplyr::select(Patient.ID,nom,PriorPARPi)%>%
  unique%>%
  ggplot(aes(x=PriorPARPi, group=as.factor(nom), fill=as.factor(nom))) +
  geom_bar(stat="count", width=0.6, position="stack")+
  xlab("Prior PARPi treatment")+
  ylab("No. of patients") +
  my_theme()+
  scale_fill_manual(values=scales::seq_gradient_pal(high = "#E64B35FF", low = "#4DBBD5FF",space="Lab")(seq(0,1,length.out=9)), 
                    name="No. of mutations") -> p.nom_parpi
p.nom_parpi

png("output/figures/nom_priorparpi.png",width=4, height=4,units="in",res=500,type="cairo")
p.nom_parpi
dev.off()

## number of mutations as a function of prior lines
df.clin %>% 
  mutate(No_PreviousLines = factor(ifelse(Number_PreviousLines >= 3, ">2",Number_PreviousLines),levels=c("1","2",">2")))%>%
  dplyr::select(Patient.ID,nom,No_PreviousLines)%>%
  unique%>%
  ggplot(aes(x=No_PreviousLines, group=as.factor(nom), fill=as.factor(nom))) +
  geom_bar(stat="count", width=0.6, position="stack")+
  xlab("No. of previous therapy lines")+
  ylab("No. of patients") +
  my_theme()+
  scale_fill_manual(values=scales::seq_gradient_pal(high = "#E64B35FF", low = "#4DBBD5FF",space="Lab")(seq(0,1,length.out=9)), 
                    name="No. of mutations") -> p.nom_prevlines
p.nom_prevlines

png("output/figures/nom_prevlines.png",width=5, height=4,units="in",res=500,type="cairo")
p.nom_prevlines
dev.off()

####compare median VAFS as a function of prior PARPi #########################
df.clin %>%
  mutate(no_lines = cut(Number_PreviousLines,breaks=c(-Inf,0,1,2,3,Inf)))%>%
  filter(CH==1)%>%
  ggboxplot(., 
            x = "no_lines",
            y = "maxVAF_CH",
            combine = TRUE,
            xlab = "Prior PARPi treatment",
            ylab = "Variant allele frequency",
            title = "",
            width = 0.3,
            #ylim = c(-4,16),
            size=0.8,
            alpha=1,
            repel=TRUE,
            #yscale = "log10",
            scales = "free",
            add = c("jitter")
  )+
  stat_compare_means()+
  scale_color_npg()+
  theme_minimal() + 
  theme(axis.title.x = element_blank()) ->p.VAF_boxplot


#####Duration of PARPi Treatment and No_previous lines vs clonesize (discrete)

my_vars_baseline=c("Number_PreviousLines",
                   "Number_PreviousPlatinumLines",
                   "no_prev_lines_binom",
                   "No_Platinum_lines_binom",
                   "Type_PreviousTherapy",
                   "PriorPARPi",
                   "Duration_PriorPARPi",
                   "Duration_PARPi_level",
                   "Age_TreatmentStartEUDARIO")

cat_vars_baseline=c("PriorPARPi",
                    "Duration_PARPi_level",
                    "No_Platinum_lines_binom")

cont_vars_baseline = setdiff(my_vars_baseline,cat_vars_baseline)

df.clin %>% CreateTableOne(strata = "CH_category",
                           vars=c(my_vars_baseline),
                           factorVars = cat_vars_baseline,
                           includeNA=FALSE,
                           #addOverall = TRUE,
                           data=.) %>% 
  print(., 
        nonnormal=cont_vars_baseline,
        exact=cat_vars_baseline,
        missing=TRUE,
        showAllLevels=TRUE,
        quote=FALSE)

##Clonesize vs no prev lines
df.clin %>%
  ggplot(.) +
  geom_mosaic(aes(x = product(CH_category, no_prev_lines_binom), fill=CH_category),alpha=1) + 
  my_theme()+
  xlab("No. of previous lines")+
  ylab("CH clone size")+
  theme(legend.position = "none")+
  scale_fill_npg()-> p.mosaic1
p.mosaic1

png("output/figures/p.mosaic1.png",width=2.5, height=2.5,units="in",res=500,type="cairo")
p.mosaic1
dev.off()


##Number of mutations vs previous lines
df.clin %>%
  mutate(nom_CH2 = cut(nom_CH,breaks=c(-Inf,0,1,2,Inf),labels=c("0","1","2",">2")))%>%
  ggplot(.) +
  geom_mosaic(aes(x = product(nom_CH2, no_prev_lines_binom), fill=nom_CH2),alpha=1) + 
  my_theme()+
  xlab("No. of previous lines")+
  ylab("No. of mutations")+
  theme(legend.position = "none")+
  scale_fill_npg()-> p.mosaic2
p.mosaic2

png("output/figures/p.mosaic2.png",width=2.5, height=2.5,units="in",res=500,type="cairo")
p.mosaic2
dev.off()


##Clonesize vs PARPi
df.clin %>%
  ggplot(.) +
  geom_mosaic(aes(x = product(CH_category, PriorPARPi), fill=CH_category),alpha=1) + 
  my_theme()+
  xlab("Prior PARPi therapy")+
  ylab("CH clone size")+
  theme(legend.position = "none")+
  scale_fill_npg()-> p.mosaic3
p.mosaic3

png("output/figures/p.mosaic3.png",width=2.5, height=2.5,units="in",res=500,type="cairo")
p.mosaic3
dev.off()

##No. mutations vs PARPi
df.clin %>%
  mutate(nom_CH2 = cut(nom_CH,breaks=c(-Inf,0,1,2,Inf),labels=c("0","1","2",">2")))%>%
  ggplot(.) +
  geom_mosaic(aes(x = product(nom_CH2, PriorPARPi), fill=nom_CH2),alpha=1) + 
  my_theme()+
  xlab("Prior PARPi therapy")+
  ylab("No. of mutations")+
  theme(legend.position = "none")+
  scale_fill_npg()-> p.mosaic4
p.mosaic4

png("output/figures/p.mosaic4.png",width=2.5, height=2.5,units="in",res=500,type="cairo")
p.mosaic4
dev.off()

png("output/figures/p.mosaic.png",width=6, height=6,units="in",res=500,type="cairo")
ggarrange(p.mosaic1,p.mosaic2,p.mosaic3,p.mosaic4)
dev.off()

####### logistic regression ########################

##DTA mutations
log.reg <- glm(DTA ~ age_dec + no_prev_lines_binom +Duration_PriorPARPi, family="binomial",data=df.clin %>% mutate(age_dec=Age_TreatmentStartEUDARIO/10))
summary(log.reg)

##odds ratios + 95%CI
data.frame(OR=exp(summary(log.reg)$coefficients[,1]),
           lCI=exp(summary(log.reg)$coefficients[,1]+ qnorm(0.025) * summary(log.reg)$coefficients[,2]),
           uCI=exp(summary(log.reg)$coefficients[,1]+ qnorm(0.975) * summary(log.reg)$coefficients[,2]))

##DDR mutations
log.reg <- glm(DDR ~ age_dec + no_prev_lines_binom+ Duration_PriorPARPi, family="binomial",data=df.clin %>% mutate(age_dec=Age_TreatmentStartEUDARIO/10, DDR = ifelse(TP53==1|PPM1D==1|ATM==1|CHEK2==1,1,0)))
summary(log.reg)

##odds ratios + 95%CI
data.frame(OR=exp(summary(log.reg)$coefficients[,1]),
           lCI=exp(summary(log.reg)$coefficients[,1]+ qnorm(0.025) * summary(log.reg)$coefficients[,2]),
           uCI=exp(summary(log.reg)$coefficients[,1]+ qnorm(0.975) * summary(log.reg)$coefficients[,2]))


