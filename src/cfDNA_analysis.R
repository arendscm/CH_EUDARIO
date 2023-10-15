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
library(ggpubr)
library(g3viz)
library(ggsci)

#####  Data preparation ####
########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')
#load('data/interim/seqdata_filtered.RData')
#load('data/interim/seqdata_filtered_cf.RData')

######## Get Patient ids
source("src/material_table.R")

######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")

##interesting gene groups
source("src/genegroup_definitions.R")

##clinical data
load("data/interim/clin.RData")

##filtered mutation data
load('data/interim/seqdata_filtered.RData')

##df preparation
df %>% mutate(gene_group =  ifelse(is.element(Gene,tp53_genes),"TP53",
                                   ifelse(is.element(Gene,hrd_genes),"HRD",
                                                 ifelse(is.element(Gene,typical_ch_genes),"CH","other myeloid")))) %>%
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

#df.filtered.c1d1 %>% 
#  filter(Patient.ID %in% (df.material %>% filter(c1d1_cf==1&c1d1_wb==1) %>% .$Patient.ID))%>%
#  filter(tag == "true",TVAF >=0.01)
#  dplyr::select(all_of(variables)) %>%
#  mutate(cfID=paste(Patient.ID,position,sep="_")) -> df.cf_wb


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
            conf.int = TRUE, 
            cor.coef = TRUE, 
            cor.method = "pearson",
            size=1,
            xlab = "cfDNA VAF", 
            ylab = "whole-blood VAF")+
  scale_y_log10(limits=c(0.001,0.5))+
  scale_x_log10(limits=c(0.001,0.5))+
  my_theme()->p.cfDNACor
p.cfDNACor

png("output/figures/p.cf.wb.png",width=4, height=4,units="in",res=500,type="cairo")
p.cfDNACor
dev.off()


###Correlation plot for all mutations found in cfDNA and WB
full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0.001,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0.001,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  filter(TVAF.y+TVAF.x>=0.008)%>%
  filter(TVAF.y > 0.2*TVAF.x,TVAF.y>0.001)%>%
  mutate(tag = ifelse(is.na(tag.x),
                      ifelse(is.na(tag.y),
                             "not tagged",
                             as.character(tag.y)),
                      as.character(tag.x)))%>%
  mutate(gene_group = ifelse(is.na(gene_group.x),gene_group.y,gene_group.x))%>%
  ggscatter(., 
            y = "TVAF.y", 
            x = "TVAF.x", 
            add = "reg.line", 
            conf.int = TRUE, 
            cor.coef = TRUE, 
            cor.method = "pearson",
            size=1,
            xlab = "cfDNA VAF", 
            ylab = "whole-blood VAF")+
  scale_y_log10(limits=c(0.001,0.5))+
  scale_x_log10(limits=c(0.001,0.5))->p.cfDNACor
p.cfDNACor

#####  Plot that shows VAF WB vs VAF ctDNA including color for group of mutation####
full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0.001,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0.001,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  filter(TVAF.y+TVAF.x>=0.008)%>%
  mutate(tag = ifelse(is.na(tag.x),
                      ifelse(is.na(tag.y),
                             "not tagged",
                             as.character(tag.y)),
                      as.character(tag.x)))%>%
  mutate(gene_group = ifelse(is.na(gene_group.x),gene_group.y,gene_group.x))%>%
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  ggplot(aes(x=TVAF.x,y=TVAF.y,
             color=gene_group
             #shape=ExonicFunc.x
             )) +
  geom_point(size=2)+
  geom_abline(slope=1,size=1,linetype=2,alpha=0.5)+
  #geom_abline(slope=1/3,size=1,linetype=3,alpha=0.5,intercept=-2)+
  scale_y_log10()+
  scale_x_log10()+
  scale_color_npg(name="Gene group",labels=c("CH genes","HR-related genes","other myeloid genes","TP53"))+
  ylab("whole-blood VAF")+
  xlab("cfDNA VAF")+
  my_theme2() -> p.cf.corr
p.cf.corr

png("output/figures/p.cf.corr.all.png",width=6.5, height=5,units="in",res=500,type="cairo")
p.cf.corr
dev.off()

### Mutationsspectrum of mutations of other than hematopoietic origin

full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  filter(TVAF.y+TVAF.x>=0.008)%>%
  mutate(tag = ifelse(is.na(tag.x),
                      ifelse(is.na(tag.y),
                             "not tagged",
                             as.character(tag.y)),
                      as.character(tag.x)))%>%
  mutate(gene_group = ifelse(is.na(gene_group.x),gene_group.y,gene_group.x))%>%
  mutate(gene_group = ifelse(gene_group == "HRD","HR-related",gene_group))%>%
  mutate(Gene = ifelse(is.na(Gene.x),Gene.y,Gene.x))%>%
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  filter(compartment == "cf")%>%
  dplyr::select(gene_group) %>% 
  table %>% 
  data.frame %>%
  ggplot(aes(x=reorder(gene_group, Freq), y=Freq, fill=gene_group)) +
  geom_bar(stat="identity", width=0.6, position = "stack")+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(position = "right")+
  ylab("No. of mutations") +
  scale_fill_npg(name="Origin",labels=c("wb" = "hematopoietic","cf" = "other"))+
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35),
        axis.ticks.y = element_blank())+
  coord_flip() -> p.mutprev_cf
p.mutprev_cf


png("output/figures/p.cf_only.mutprev.png",width=5, height=3,units="in",res=500,type="cairo")
p.mutprev_cf
dev.off()

#####   Mutationspectrum for mutations in WB and cf only mutations by group------------------
full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  filter(TVAF.y+TVAF.x>=0.008)%>%
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
  ylab("No. of mutations") +
  scale_fill_npg(name="Origin",labels=c("cf" = "other", "wb" = "hematopoietic"))+
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35),
        axis.ticks.y = element_blank())+
  coord_flip() -> p.mutprev
p.mutprev


png("output/figures/p.cf.mutprev.png",width=5, height=3,units="in",res=500,type="cairo")
p.mutprev
dev.off()

##numbers

full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  filter(TVAF.y+TVAF.x>=0.008)%>%
  mutate(tag = ifelse(is.na(tag.x),
                      ifelse(is.na(tag.y),
                             "not tagged",
                             as.character(tag.y)),
                      as.character(tag.x)))%>%
  mutate(gene_group = ifelse(is.na(gene_group.x),gene_group.y,gene_group.x))%>%
  mutate(Gene = ifelse(is.na(Gene.x),Gene.y,Gene.x))%>%
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  dplyr::select(gene_group, compartment) %>% table

### Mutationspectrum cf vs wb by gene 
full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  filter(TVAF.y+TVAF.x>=0.008)%>%
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
  group_by(Gene) %>%
  mutate(Freq_all =sum(Freq))%>%
  data.frame%>%
  filter(Freq_all > 1)%>%
  ggplot(aes(x=reorder(Gene, Freq_all), y=Freq, fill=compartment)) +
  geom_bar(stat="identity", width=0.6, position = "stack")+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(position = "right")+
  ylab("No. of mutations") +
  scale_fill_npg(name="Origin",labels=c("wb" = "hematopoietic","cf" = "other"),breaks=c("wb","cf"))+
  my_theme() +
  theme(axis.text.y=element_text(angle=0,
                                 hjust=1,
                                 vjust=0.35,
                                 face="italic"),
        axis.ticks.y = element_blank(),
        legend.position = c(0.8,0.2))+
  coord_flip() -> p.mutprev_all
p.mutprev_all


png("output/figures/p.cf.mutprev_all.png",width=4, height=6,units="in",res=500,type="cairo")
p.mutprev_all
dev.off()

################## Correlation of cfDNA VAF with tumorburden and CA125
library(GGally)
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
  mutate(compartment = ifelse(TVAF.x > TVAF.y*10,"cf","wb"))%>%
  filter(gene_group=="TP53")%>%
  filter(compartment=="cf")%>%
  group_by(Patient.ID.x)%>%
  mutate(maxVAF = max(TVAF.x))%>%
  data.frame%>%
  filter(TVAF.x == maxVAF)%>%
  mutate(Patient.ID = Patient.ID.x,
         TP53_mut = AAChange.x,
         VAF_TP53 = log10(TVAF.x))%>%
  dplyr::select(Patient.ID,TP53_mut,VAF_TP53)%>%
  left_join(.,df.clin,by=c("Patient.ID"))%>%
  mutate(logCA=log10(CA125))%>%
  dplyr::select(logCA,TumorBurden_baseline,VAF_TP53)%>%
  ggpairs()


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
         Patient.ID = Patient.ID.x,
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
  dplyr::select(Gene,AAChange,Patient.ID,Func,ExonicFunc,Chr,Start,End,Ref,Alt,TVAF,tag,COSMIC)-> df.tp53

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

##Save RData for further use
full_join(df.cf%>% filter(Patmut %in% Patmut_all),df.cf_wb %>% filter(Patmut %in% Patmut_all),by="cfID") %>% 
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0.001,TVAF.y))%>%
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0.001,TVAF.x))%>%
  filter(TR2.x>9|TR2.y>9)%>%
  filter(TVAF.y+TVAF.x>=0.008)%>%
  mutate(tag = ifelse(is.na(tag.x),
                      ifelse(is.na(tag.y),
                             "not tagged",
                             as.character(tag.y)),
                      as.character(tag.x)))%>%
  mutate(gene_group = ifelse(is.na(gene_group.x),gene_group.y,gene_group.x))%>%
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb")) -> df.cf_wb_c1d1

tempdata <-ls()
rm(list=tempdata[!is.element(tempdata,c("df.cf_wb_c1d1"))])
rm(tempdata)

save.image("data/interim/cf_wb.RData")




