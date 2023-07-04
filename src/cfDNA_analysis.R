# ______________________________________________________________________________
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: cfDNA Mutation analysis 
#
# Input: df from data/interim/seqdata.RData
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
library(g3viz)
library(ggsci)

#####  Data preparation ####
########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')
load('data/interim/seqdata_filtered.RData')
load('data/interim/seqdata_filtered_cf.RData')

######## Get Patient ids
source("src/material_table.R")

######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")

##interesting gene groups
source("src/genegroup_definitions.R")

##df preparation
df %>% mutate(gene_group =  ifelse(is.element(Gene,tp53_genes),"TP53",
                                   ifelse(is.element(Gene,ppm1d_genes),"PPM1D",
                                          ifelse(is.element(Gene,hrd_genes),"HRD",
                                                 ifelse(is.element(Gene,ch_genes_without_HRD),"CH","other"))))) %>%
  mutate(frac_mut = mutFreq/n.lane) -> df

##relevant variables
variables <- c("Patient.ID","Sample_orig","mutID","position","Sample", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "GeneDetail", "ExonicFunc", "AAChange", "cytoBand","readDepth", "TR1", "TR1_plus", "TR1_minus", "TR2", "TR2_plus", "TR2_minus", "TVAF", "AF", "avsnp150","cosmic92_coding","snp","mutFreq","p.binom","n.mut","n.material","sum_cf","sum_wb","Material","tag", "Patmut","gene_group","frac_mut")


#dataframe with all mutation calls from cfDNA at c1d1          
df %>% 
  filter(is.na(replicate))%>%
  filter(Patient.ID %in% (df.material %>% filter(c1d1_cf==1&c1d1_wb==1) %>% .$Patient.ID))%>%
  filter(Material=="cf") %>% 
  filter(Visite == "C1D1")%>%
  dplyr::select(all_of(variables)) %>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))-> df.cf

#data frame with all mutation calls from WB samples that have matched cfDNA samples
df %>% 
  filter(is.na(replicate))%>%
  filter(Patient.ID %in% (df.material %>% filter(c1d1_cf==1&c1d1_wb==1) %>% .$Patient.ID))%>%
  filter(Material=="wb") %>% 
  filter(Visite == "C1D1")%>%
  dplyr::select(all_of(variables)) %>%
  mutate(cfID=paste(Patient.ID,position,sep="_")) -> df.cf_wb

#Mutations to work with
df.cf %>% 
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.01)%>%
  filter(snp==FALSE)%>%
  filter(tag!="germline")%>%
  filter(tag!="false")%>%
  filter(TVAF >= 0.008)%>%
  filter(p.binom < -12)%>%
  filter(frac_mut < 0.2)%>%
  full_join(.,df.cf %>% filter(tag == "true"|tag=="cf-only")) %>% .$Patmut -> Patmut_cf

df.cf_wb %>% 
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.01)%>%
  filter(snp==FALSE)%>%
  filter(tag!="germline")%>%
  filter(tag!="false")%>%
  filter(TVAF >= 0.008)%>%
  filter(p.binom < -12)%>%
  filter(frac_mut < 0.2)%>%
  full_join(.,df.cf_wb %>% filter(tag == "true")) %>% .$Patmut -> Patmut_wb

Patmut_all <- c(Patmut_cf,Patmut_wb)%>% unique

#####  Identity check via SNP  ####
left_join(df.cf,df.cf_wb,by="cfID")%>%
  filter(snp.x == 1) %>% 
  ggplot(aes(x=Patient.ID.x,y=TVAF.x-TVAF.y)) +
  geom_point()+
  coord_flip()+
  theme(axis.text.x =element_text(angle=90)) ->p.cf.snp

p.cf.snp

png("output/figures/p.cf.snp.png",width=5, height=7,units="in",res=500,type="cairo")
p.cf.snp
dev.off()

#####  Correlation Plot WB vs cfDNA calls tagged true, only calls that are present in WB ####
df.cf_wb %>%   
  filter(tag == "true")%>%
  filter(TVAF >= 0.01)%>%
  left_join(.,df.cf,by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0.001,TVAF.y))%>%
  ggscatter(., 
            y = "TVAF.x", 
            x = "TVAF.y", 
            add = "reg.line", 
            #conf.int = TRUE, 
            cor.coef = TRUE, 
            cor.method = "pearson",
            size=1,
            xlab = "cfDNA VAF", 
            ylab = "whole-blood VAF")+
  scale_y_log10()+
  scale_x_log10()->p.cfDNACor
p.cfDNACor

png("output/figures/p.cf.wb.png",width=4, height=4,units="in",res=500,type="cairo")
p.cfDNACor
dev.off()


#####  Plot that shows VAF WB vs VAF ctDNA including color for group of mutation - NO Filter####
full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0.001,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0.001,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  mutate(tag = ifelse(is.na(tag.x),
                      ifelse(is.na(tag.y),
                             "not tagged",
                             as.character(tag.y)),
                      as.character(tag.x)))%>%
  mutate(gene_group = ifelse(is.na(gene_group.x),gene_group.y,gene_group.x))%>%
  ggplot(aes(x=TVAF.x,y=TVAF.y,
             color=gene_group,
             #shape=ExonicFunc.x
             )) +
  geom_point(size=2)+
  geom_abline(slope=1,size=1,linetype=2,alpha=0.5)+
  geom_abline(slope=1/3,size=1,linetype=3,alpha=0.5,intercept=-2)+
  scale_y_log10()+
  scale_x_log10()+
  scale_color_npg()+
  ylab("whole blood VAF")+
  xlab("cfDNA VAF")+
  theme_minimal()->p.cf.corr
p.cf.corr

png("output/figures/p.cf.corr.all.png",width=10, height=6,units="in",res=500,type="cairo")
p.cf.corr
dev.off()

#####  Max: Mutationspectrum for mutations in WB and cf only mutations by group------------------
full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  mutate(tag = ifelse(is.na(tag.x),
                      ifelse(is.na(tag.y),
                             "not tagged",
                             as.character(tag.y)),
                      as.character(tag.x)))%>%
  mutate(gene_group = ifelse(is.na(gene_group.x),gene_group.y,gene_group.x))%>%
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  dplyr::select(gene_group, compartment) %>% 
  table %>% 
  data.frame %>%
  ggplot(aes(x=reorder(gene_group, Freq), y=Freq, fill=compartment)) +
  geom_bar(stat="identity", width=0.6, position = "stack")+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(position = "right")+
  ylab("Gene Mutation Frequency") +
  scale_fill_npg()+
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35),
        axis.ticks.y = element_blank())+
  coord_flip() -> p.mutprev
p.mutprev

### Mutationspectrum cf vs wb by gene 
full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  mutate(tag = ifelse(is.na(tag.x),
                      ifelse(is.na(tag.y),
                             "not tagged",
                             as.character(tag.y)),
                      as.character(tag.x)))%>%
  mutate(gene_group = ifelse(is.na(gene_group.x),gene_group.y,gene_group.x))%>%
  mutate(Gene = ifelse(is.na(Gene.x),Gene.y,Gene.x))%>%
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  dplyr::select(Gene, compartment) %>% 
  table %>% 
  data.frame  %>%
  ggplot(aes(x=reorder(Gene, Freq), y=Freq, fill=compartment)) +
  geom_bar(stat="identity", width=0.6, position = "stack")+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(position = "right")+
  ylab("Gene Mutation Frequency") +
  scale_fill_npg()+
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35),
        axis.ticks.y = element_blank())+
  coord_flip() -> p.mutprev
p.mutprev


##Lolliplots for TP53 mutations
full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  mutate(tag = ifelse(is.na(tag.x),
                      ifelse(is.na(tag.y),
                             "not tagged",
                             as.character(tag.y)),
                      as.character(tag.x)))%>%
  mutate(gene_group = ifelse(is.na(gene_group.x),gene_group.y,gene_group.x))%>%
  mutate(Gene = ifelse(is.na(Gene.x),Gene.y,Gene.x))%>%
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  filter(Gene=="TP53")%>%
  mutate(AAChange = AAChange.x,
         Sample = Patient.ID.x,
         ExonicFunc=ExonicFunc.x,
         Chr=Chr.x,
         Start=Start.x,
         End=End.x,
         Ref=Ref.x,
         Alt=Alt.x,
         TVAF=TVAF.x,
         tag = tag.x,
         Func=Func.x,
         COSMIC=cosmic92_coding.x)%>%
  dplyr::select(Gene,AAChange,Sample,Func,ExonicFunc,Chr,Start,End,Ref,Alt,TVAF,tag,COSMIC)-> df.tp53

df.tp53_wb <- df.tp53 %>% filter(tag=="true")%>%makeMAF
df.tp53_cf <- df.tp53 %>% filter(tag=="cf-only")%>%makeMAF
  
plot.options <- g3Lollipop.theme(theme.name = "cbioportal",
                                 title.text = "TP53 CH",
                                 y.axis.label = "# of Mutations")

g3Lollipop(df.tp53_wb,
           gene.symbol = "TP53",
           btn.style = "gray", # gray-style chart download buttons
           plot.options = plot.options,
           factor.col = "Variant_Type",
           save.png.btn	= FALSE,
           save.svg.btn = FALSE,
           output.filename = "cbioportal_theme")


plot.options <- g3Lollipop.theme(theme.name = "cbioportal",
                                 title.text = "TP53 cfDNA",
                                 y.axis.label = "# of Mutations")

g3Lollipop(df.tp53_cf,
           gene.symbol = "TP53",
           btn.style = "gray", # gray-style chart download buttons
           plot.options = plot.options,
           factor.col = "Variant_Type",
           save.png.btn	= FALSE,
           save.svg.btn = FALSE,
           output.filename = "cbioportal_theme")



