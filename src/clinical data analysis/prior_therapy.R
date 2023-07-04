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
                   "No_Platinum_lines_binom",
                   "Type_PreviousTherapy",
                   "PriorPARPi",
                   "Duration_PriorPARPi",
                   "Age_TreatmentStartEUDARIO")

cat_vars_baseline=c("Type_PreviousTherapy",
                    "PriorPARPi",
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


####### logistic regression ########################

log.reg <- glm(CH ~ Age_TreatmentStartEUDARIO  +PriorPARPi + No_Platinum_lines_binom, family="binomial",data=df.clin)
summary(log.reg)


##### compare number of mutations as a function of prior PARPi

df.clin %>% 
  dplyr::select(Patient.ID,nom,PriorPARPi)%>%
  unique%>%
  ggplot(aes(x=PriorPARPi, group=as.factor(nom), fill=as.factor(nom))) +
  geom_bar(stat="count", width=0.6, position="stack")+
  xlab("Prior PARPi")+
  ylab("Patient count") +
  my_theme()+
  scale_fill_manual(values=scales::seq_gradient_pal(high = "#E64B35FF", low = "#4DBBD5FF",space="Lab")(seq(0,1,length.out=9)), 
                    name="No. of mutations") -> p.nom_parpi
p.nom_parpi

png("output/figures/nom_priorparpi.png",width=6, height=6,units="in",res=500,type="cairo")
p.nom_parpi
dev.off()

## number of mutations as a function of prior lines
df.clin %>% 
  mutate(No_PreviousLines = factor(ifelse(Number_PreviousLines >= 3, ">2",Number_PreviousLines),levels=c("1","2",">2")))%>%
  dplyr::select(Patient.ID,nom,No_PreviousLines)%>%
  unique%>%
  ggplot(aes(x=No_PreviousLines, group=as.factor(nom), fill=as.factor(nom))) +
  geom_bar(stat="count", width=0.6, position="stack")+
  xlab("No. previous therapy lines")+
  ylab("Patient count") +
  my_theme()+
  scale_fill_manual(values=scales::seq_gradient_pal(high = "#E64B35FF", low = "#4DBBD5FF",space="Lab")(seq(0,1,length.out=9)), 
                    name="No. of mutations") -> p.nom_prevlines
p.nom_prevlines

png("output/figures/nom_prevlines.png",width=6, height=6,units="in",res=500,type="cairo")
p.nom_prevlines
dev.off()

####compare median VAFS as a function of prior PARPi #########################
df.clin %>%
  filter(CH==1)%>%
  ggboxplot(., 
            x = "Number_PreviousLines",
            y = "maxVAF",
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

####compare VAFS as a function of prior PARPi time#########################
