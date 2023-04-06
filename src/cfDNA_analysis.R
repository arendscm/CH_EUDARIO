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
###### Dependencies   #####
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
###### Data preparation ####
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
Categories<-c('CH','HRD','other', 'TP53')
#data frame with mutation calls from cfDNA             
df %>% 
  filter(Material=="cf") %>% 
  filter(Visite == "C1D1")%>%
  filter(is.na(replicate))%>%
  filter(!is.element(Sample.ID, failedSamples))%>%
  dplyr::select(all_of(variables)) %>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))-> df.cf
df %>% 
  filter(Material=="wb") %>% 
  filter(Visite == "C1D1")%>%
  filter(is.na(replicate))%>%
  filter(!is.element(Sample.ID, failedSamples))%>%
  dplyr::select(all_of(variables)) %>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))-> df.wb

#anti_join(df.cf, df.wb, by= "cfID")->df.cf_only
#save(df.cf_only,file="data/interim/df.cf_only.RDATA")

#data frame with mutation calls from WB samples that have matched cfDNA samples
df %>% 
  filter(is.element(Patient.ID,df.cf$Patient.ID)) %>% 
  filter(is.na(replicate))%>%
  filter(Material=="wb") %>% 
  filter(Visite == "C1D1")%>%
  filter(!is.element(Sample.ID, failedSamples))%>%
  dplyr::select(all_of(variables)) %>%
  mutate(cfID=paste(Patient.ID,position,sep="_")) -> df.cf_wb

#####  identity check via SNP  ####
left_join(df.cf,df.cf_wb,by="cfID")%>%
  filter(snp.x == 1) %>% 
  ggplot(aes(x=Patient.ID.x,y=TVAF.x-TVAF.y)) +
  geom_point()+coord_flip()->p.cf.snp
png("output/figures/p.cf.snp.png",width=5, height=7,units="in",res=500,type="cairo")
p.cf.snp
dev.off()


#####  Correlation Plot WB vs cfDNA all calls, only calls that are present in WB ####
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x))%>%
  filter(TVAF.y < 0.4&TVAF.x<0.4)%>%
  filter(!snp.y&!snp.x)%>%
  ggscatter(., 
            x = "TVAF.x", 
            y = "TVAF.y", 
            add = "reg.line", 
            #conf.int = TRUE, 
            cor.coef = TRUE, 
            cor.method = "pearson",
            size=1,
            ylab = "VAF wholeblood", 
            xlab = "VAF cfDNA")->p.cfDNACor
p.cfDNACor

#####  Correlation Plot WB vs cfDNA calls tagged true, only calls that are present in WB ####
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(tag.y == "true")%>%
  filter(TVAF.y > 0.01)%>%
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
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  #filter(!is.element(Sample.x,mismatch))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes),"CH",
                       ifelse(is.element(Gene.x,tp53_genes),"TP53",
                              ifelse(is.element(Gene.x,hrd_genes),"HRD",
                                     ifelse(is.element(Gene.x,brca_genes),"BRCA",
                                            ifelse(is.element(Gene.x,ppm1d_genes),"PPM1D","other"))))))%>%
  #filter(gene != "other") %>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding.x,"ovary")) %>%
  #filter(p.binom.x <= -10) %>%
  #filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  #filter(ExonicFunc.x != "synonymous SNV")%>%
  #filter(AF.x<0.1)%>%
  #filter(snp.x==FALSE)%>%
  #filter(TVAF.x>0.005|TVAF.y>0.005)%>%
  #filter(TR2.y > 19|TR2.x>19)%>%
  ggplot(aes(x=TVAF.x,y=TVAF.y,
             color=gene,
             #shape=ExonicFunc.x
             )) +
  geom_point(size=2)+
  geom_abline(slope=1,size=1,linetype=2,alpha=0.5)+
  geom_abline(slope=1/3,size=1,linetype=3,alpha=0.5,intercept=-2)+
  scale_y_log10()+
  scale_x_log10()+

  #facet_wrap(~ Patient.ID.x, ncol=4, dir="h")+
  scale_color_viridis(discrete=TRUE)+
  ylab("whole blood VAF")+
  xlab("cfDNA VAF")+
  theme_minimal()->p.cf.corr

png("output/figures/p.cf.corr.all.png",width=10, height=6,units="in",res=500,type="cairo")
p.cf.corr
dev.off()

#####  Plot that shows VAF WB vs VAF ctDNA including color for group of mutation- Filter 1 ####
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  #filter(!is.element(Sample.x,mismatch))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes),"CH",
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
  ggplot(., aes(x=TVAF.x,y=TVAF.y,
             color=gene,
             #shape=ExonicFunc.x
             )) +
  geom_point(size=2)+
  geom_abline(slope=1,size=1,linetype=2,alpha=0.5)+
  geom_abline(slope=1/3,size=1,linetype=3,alpha=0.5,intercept=-2)+
  scale_y_log10()+
  scale_x_log10()+
  #facet_wrap(~ Patient.ID.x, ncol=4, dir="h")+
  scale_color_viridis(discrete=TRUE)+
  ylab("whole-blood VAF")+
  xlab("cfDNA VAF")+
  theme_minimal()->p.cf.corr
p.cf.corr

<<<<<<< HEAD
png("output/figures/p.cf.corr.filter1.png",width=6, height=4,units="in",res=500,type="cairo")
=======
png("output/figures/p.cf.corr.filter1.png",width=10, height=10,units="in",res=500,type="cairo")
>>>>>>> a68ef340685e7eed40598bfdd34819df7b3e8e52
p.cf.corr
dev.off()


#### below line composition analysis
load("data/interim/seqdata_filtered_cf.RData")

df_filtered_cf%>%
  filter(tag == "true" | tag == "tumour")%>%
  mutate(gene = ifelse(is.element(Gene, ch_genes_without_HRD), "CH",
                       ifelse(is.element(Gene, tp53_genes), "TP53",
                              ifelse(is.element(Gene, hrd_genes), "HRD",
                                     ifelse(is.element(Gene, brca_genes), "BRCA",
                                            ifelse(is.element(Gene, ppm1d_genes), "PPM1D", "other"))))))%>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))%>%
  left_join(.,df.cf_wb, by= "cfID")%>%
  filter(TVAF.y < TVAF.x/3 | is.na(TVAF.y))%>%
  filter(TVAF.y < 0.001| is.na(TVAF.y))-> df.cf_only


# Create a table of gene counts
gene_counts <- table(df.cf_only$gene)

# Calculate the percentage of variants in each gene category
gene_percents <- prop.table(gene_counts) * 100

# Convert the table object to a data frame
gene_percents_df <- as.data.frame(gene_percents)

# Create a bar plot of the percentage of variants in each gene category
ggplot(gene_percents_df, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", fill = "#486081") +
  labs(title = "Percentage of Variants by Gene Category",
       x = "Gene Category",
       y = "Percentage")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.35, face = "italic", size = 16),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.35, face = "italic", size = 16),
        plot.title = element_text(size = 20, face = "bold")) ->p.cfonly.category

png("output/figures/p.cfonly.category.png",width=10, height=10,units="in",res=500,type="cairo")
p.cfonly.category
dev.off()


##PLOTs
nop <- ids%>%
  filter(Visite == "C1D1" & Material == "cf")%>%
  filter(!is.element(Sample.ID, failedSamples))%>%
  select(.,Patient.ID)%>%
  unique()%>%nrow

########   Gene Mutation Prevalence Plot (plots number of gene-x-mutated patients)  #####
for (x in Categories)
{
  df.cf_only%>% 
    filter(gene == x)%>%
    dplyr::select(Sample.x, Gene.x)%>%
    data.frame %>% 
    unique %>% 
    dplyr::select(Gene.x) %>% 
    table %>% 
    data.frame %>% 
    filter(Freq>0) %>% 
    mutate(prev = Freq/nop) %>% 
    arrange(prev) -> prev.table
  names(prev.table)<- c("Gene","Freq","prev")
  
  prev.table %>%
    ggplot(aes(x = reorder(Gene, Freq), y = prev)) +
    geom_bar(stat = "identity", width = 0.6, fill = "#486081") +
    geom_text(aes(label = Freq), hjust = -1, vjust = 0.35, size = 4) +
    xlab("") +
    scale_y_continuous(labels = percent, limits = c(0, 0.6), position = "right") +
    ylab("Gene Mutation Prevalence [%]") +
    my_theme() +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.35, face = "italic", size = 16),
          axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.35, face = "italic", size = 14)) +
    coord_flip()  -> p.mutprev.composition
  
  
  png(paste0("output/figures/",x,"_cfonly_composition.png"),
      width=10,
      height=10,
      units="in",
      res=500,
      type="cairo")
  print(p.mutprev.composition)
  dev.off()
}



##### Correlation plot Filter 1 - TP53 mutations only 
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(!is.element(Patient.ID.x,c("4202008","4207011","4220002","4221004","4221032")))%>%
  #filter(!is.element(Sample.x,mismatch))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes),"CH",
                       ifelse(is.element(Gene.x,tp53_genes),"TP53",
                              ifelse(is.element(Gene.x,hrd_genes),"HRD",
                                     ifelse(is.element(Gene.x,brca_genes),"BRCA",
                                            ifelse(is.element(Gene.x,ppm1d_genes),"PPM1D","other"))))))%>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding.x,"ovary")) %>%
  filter(p.binom.x == -Inf|p.binom.y == -Inf) %>%
  filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  filter(ExonicFunc.x != "synonymous SNV")%>%
  filter(AF.x<0.1)%>%
  #filter(snp.x==FALSE)%>%
  filter(TR2.x>19|TR2.y>19) %>% 
  filter(TVAF.x>0.001|TVAF.y>0.001)%>%
  filter(gene=="TP53")%>%
  #filter(cosmic_ovary)%>%
  ggplot(aes(x=TVAF.x,y=TVAF.y,color=cosmic_ovary))+
  geom_point(size=4)+
  geom_abline(slope=1)+
  scale_color_viridis(discrete=TRUE)+
  scale_x_log10(limits=c(0.0005,0.5)) +
  scale_y_log10(limits=c(0.0005,0.5)) +
  theme_minimal() -> p.TP53_cosmic_cf_wb

png("output/figures/p.TP53_cosmic_cf_wb.png",width=10, height=6,units="in",res=500,type="cairo")
p.TP53_cosmic_cf_wb
dev.off()

#####  detecting BRCA mutations in cfDNA (this will later on also be important when looking for BRCA reversion mutations)####
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  #filter(!is.element(Sample.x,mismatch))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes),"CH",
                       ifelse(is.element(Gene.x,tp53_genes),"TP53",
                              ifelse(is.element(Gene.x,hrd_genes),"HRD",
                                     ifelse(is.element(Gene.x,brca_genes),"BRCA",
                                            ifelse(is.element(Gene.x,ppm1d_genes),"PPM1D","other"))))))%>%
  filter(gene != "other") %>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding.x,"ovary")) %>%
  filter(p.binom.x < -10) %>%
  filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  filter(ExonicFunc.x != "synonymous SNV")%>%
  filter(AF.x<0.1)%>%
  #filter(snp.x==FALSE)%>%
  filter(TR2.x>9|TR2.y>9) %>% 
  filter(TVAF.x>0.005|TVAF.y>0.005)%>%
  filter(gene=="BRCA")%>%
  ggplot(aes(x=TVAF.x,y=TVAF.y,color=gene))+
  geom_point(size=4)+
  geom_abline(slope=1)+
  scale_color_viridis(discrete=TRUE)+
  scale_x_log10(limits=c(0.0005,0.5)) +
  scale_y_log10(limits=c(0.0005,0.5)) +
  theme_minimal()

#####  Determine Overlap CH and HRD in plasma samples ####
df.filtered%>%
  filter(Material=="cf") %>% 
  filter(Visite == "C1D1")%>%
  filter(is.na(replicate))%>%
  filter(TVAF >= 0.01 & TVAF <= 0.3)->df.filtered.cf

nop<-df.filtered.cf%>%
  select(.,Patient.ID)%>%
  unique()
source("src/global_functions_themes.R")

df_filtered_cf %>% 
  #filter(Gene %in% ch_genes)%>%  #only CH panel, when we say: this is the prevalence plot for CH in these patients
  dplyr::select(Sample, Gene) %>% 
  data.frame %>% 
  unique %>% 
  dplyr::select(Gene) %>% 
  table %>% 
  data.frame %>% 
  filter(Freq >0) %>% 
  mutate(prev = Freq/nop) %>% 
  arrange(prev) -> prev.table
names(prev.table)<- c( "Gene","Freq","prev")

prev.table  %>%
  mutate(HRD = ifelse(is.element(Gene,hrd_genes),"HRD","non HRD"))%>%
  ggplot(aes(x=reorder(Gene, Freq), y=prev, fill=HRD)) +
  geom_bar(stat="identity", width=0.6)+
  geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(labels = percent,limits=c(0,0.35), position = "right")+
  ylab("Gene Mutation Prevalence [%]") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  coord_flip() + 
  scale_fill_manual(values = c("non HRD" = "#486081", "HRD" = "#88acd4")) -> p.mutprev
p.mutprev

###### serial cf data exploration 2 timepoints cf only ####
df %>% 
  filter(Material=="cf")%>%
  filter(c1d1_cf==1&eot_cf==1)%>%  
  #filter(c1d1_cf==1&eot_cf==1 & c7d1_cf == 0)%>% #for subsetting patients with 2 or 3 timepoints
  filter(n.material>=2)%>% 
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.01)%>%
  filter(snp==FALSE)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA","other")))))%>%
  filter(maxVAF > 0.005,
         maxTR2 > 7,
         minpbinom < -8) %>%
  #filter(TVAF < 0.38) %>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=position,color=gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=6, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.cf.serial

png("output/figures/p.cf.serial.png",width=10, height=5,units="in",res=500,type="cairo")
p.cf.serial
dev.off()

#####  preparation: find mutations that are present in both wb and cf  #####
df%>%  mutate(cfID=paste(Patient.ID,position,sep="_"))%>%
  filter(Material=="wb")%>%
  group_by(Visite)%>%
  dplyr::select(cfID)->wb.cfID

df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(Material=="cf")%>%
  group_by(Visite)%>%
  mutate(wb=is.element(mutID,wb.mutID))->cf.cfID

#####  serial cf data exploration 2 timepoints cf and wb - wb=full line ; plasma=dashed ####
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=gene,group=pat_material,shape=Material),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,color=gene,linetype=Material),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=6, scales="free", dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial


#####  facet by patient and material including "cf_only" ####
##serial cf data exploration 2 timepoints cf and wb
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n()) %>% data.frame()%>% filter(n.visit==1)%>% filter(Material=="cf") %>% mutate(Material="cf_only") ->df.cfonly

df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  full_join(.,df.cfonly)%>%
  filter(is.element(Patient.ID,c("4202011","4207004","4202001","4221010")))%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  filter(Patient.ID %in% c("4203009","4203002", "4221004", "4221006", "4221023"))%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=gene,group=pat_material),size=1,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,color=gene),size=0.5,na.rm=FALSE) + 
  facet_grid(Material ~ Patient.ID) +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial

png("output/figures/p.cf.serial_cf,cf_only,wb_Patientsexample.png",width=6, height=4,units="in",res=500,type="cairo")
p.cf.serial
dev.off()


#####  Test single patient by mutation ####
df %>% 
  filter(Patient.ID=="4221010")%>%
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  mutate(AAChange_mod = str_replace(AAChange, ":NM.*:p.", ":p."))%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=Gene,group=pat_material,shape=Material),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,color=Gene,linetype=Material),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ AAChange_mod, ncol=6, dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  #scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial

##serial cf data exploration 2 timepoints cf and wb TP53 only
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(Gene=="TP53")%>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding,"ovary")) %>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  filter(maxVAF > 0.001,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=factor(n.visit),group=pat_material,shape=Material),size=2,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,linetype=Material),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=5, dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="detected in WB/cfDNA") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial


##### serial cf data exploration 2 timepoints cf and wb HRD only
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(is.element(Gene,hrd_genes))%>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding,"ovary")) %>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  filter(maxVAF > 0.005,
         maxTR2 > 19,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=factor(n.visit),group=pat_material,shape=Material),size=2,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,linetype=Material,color=cosmic_ovary),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=5, dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="detected in WB/cfDNA") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()-> p.cf.serial


###test cf serial with BRCA mutations
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding,"ovary")) %>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  filter(gene=="BRCA")%>%
  filter(maxVAF > 0.005,
         maxTR2 > 19,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=factor(n.visit),group=pat_material,shape=Material),size=2,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=pat_material_pos,linetype=Material,color=cosmic_ovary),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=5, dir="h") +
  #scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="detected in WB/cfDNA") +
  scale_y_log10(limits=c(0.0005,0.5))+
  theme_minimal()

#####  just for fun: dynamics cfVAF vs wbVAF ####
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  filter(Material == "wb")%>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))->df1
df %>% 
  filter(c1d1_cf==1&eot_cf==1)%>%
  filter(c1d1_wb==1&eot_wb==1)%>%
  filter(is.na(replicate))%>%
  mutate(cfID_visit = paste(Patient.ID,position,Visite,sep="_"))%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA",
                                            ifelse(is.element(Gene,ppm1d_genes),"PPM1D","other"))))))%>%
  group_by(cfID_visit)%>%
  mutate(n.visit=n())%>%
  data.frame%>%
  filter(n.material>=2|n.visit>1)%>%
  filter(AF<0.01)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         minVAF = min(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(pat_material = paste(Patient.ID,Material,sep="_"))%>%
  mutate(pat_material_pos = paste(Patient.ID,Material,position,sep="_"))%>%
  filter(maxVAF > 0.01,
         maxTR2 > 9,
         minpbinom < -10) %>%
  filter(maxVAF<0.4|minVAF<0.4) %>%
  filter(TVAF<0.35)%>%
  filter(Material == "cf")%>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))->df2

full_join(df1,df2,by=c("cfID_visit"))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  filter(TVAF.x!=0)%>%
  ggplot(aes(x=log(TVAF.x+0.001),y=log10(TVAF.y+0.001),
             color=gene.x)) +
  geom_point(aes(shape=Visite.y),size=3)+
  geom_line(aes(group=Patmut.y),size=1)+
  #scale_y_log10()+
  #scale_x_log10()+
  scale_color_viridis(discrete=TRUE)+
  ylab("whole blood VAF")+
  xlab("cfDNA VAF")+
  theme_minimal()



#png("output/figures/p.cf.serial.png",width=6, height=4,units="in",res=500,type="cairo")
#p.cf.serial
#dev.off()


df %>% 
  filter(Material=="cf")%>%
  filter(c1d1_cf==1&c7d1_cf==1&eot_cf==1)%>%
  filter(Patient.ID %in% c("4203001", "4203004"))%>%
  filter(n.material==3)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.01)%>%
  filter(snp==FALSE)%>%
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF),
         maxTR2 = max(TR2),
         minpbinom = min(p.binom)) %>%
  data.frame()%>%
  mutate(gene = ifelse(is.element(Gene,ch_genes),"CH",
                       ifelse(is.element(Gene,tp53_genes),"TP53",
                              ifelse(is.element(Gene,hrd_genes),"HRD",
                                     ifelse(is.element(Gene,brca_genes),"BRCA","other")))))%>%
  filter(maxVAF > 0.005,
         maxTR2 > 7,
         minpbinom < -8)%>%
  .$Patmut->Patmut.serial

df%>%
  filter(Patmut %in% Patmut.serial)%>%
  filter(Material == "cf")%>%
  #filter(TVAF < 0.38) %>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=Gene,group=Patient.ID),size=1,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=position,color=Gene),size=0.5,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=2, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.3)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal() -> p.cf.serial

png("output/figures/p.cf.serialC7.png",width=6, height=4,units="in",res=500,type="cairo")
p.cf.serial
dev.off()

#####  Lolliplot for TP53 muts------------------------------------------------------------------------
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(!is.element(Sample.x,mismatch))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes),"CH",
                       ifelse(is.element(Gene.x,tp53_genes),"TP53",
                              ifelse(is.element(Gene.x,hrd_genes),"HRD",
                                     ifelse(is.element(Gene.x,brca_genes),"BRCA","other")))))%>%
  filter(gene != "other") %>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding.x,"ovary")) %>%
  filter(mutFreq.x < 10) %>%
  filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  filter(ExonicFunc.x != "synonymous SNV")%>%
  filter(AF.x<0.1)%>%
  #filter(snp.x==FALSE)%>%
  filter(TR2.x>9|TR2.y>9) %>% 
  filter(TVAF.x>0.005|TVAF.y>0.005)%>%
  filter(gene=="TP53") %>% 
  mutate(origin = ifelse(TVAF.y>0,"WB","ctDNA"))%>%
  filter(Func.x == "exonic")%>%
  separate(.,AAChange.x,
           into=c("transcript1","rest"),
           sep=",",
           remove=TRUE,
           convert=FALSE)%>%
  separate(.,transcript1,
           into=c("gene","transcript_name","exon","DNAchange","amino_acid_change"),
           sep = ":",
           remove = TRUE,
           convert = FALSE)%>%
  #mutate(ExonicFunc = ifelse(ExonicFunc=="nonsynonymous SNV","Missense","Truncating"))%>%
  mutate(AA_pos.x = as.character(amino_acid_change))%>%
  separate(.,AA_pos.x,
           into=c("p","AA_pos1"),
           sep = "p.",
           remove = TRUE,
           convert = FALSE)%>%
  separate(.,AA_pos1,
           into=c("AA_pos2","rest"),
           sep = "fs",
           remove = TRUE,
           convert = FALSE)%>%
  mutate(AA_position = extract_numeric(AA_pos2))%>%
  dplyr::select(gene,amino_acid_change,AA_position,Sample.x,ExonicFunc.x,Chr.x,Start.x,End.x,Ref.x,Alt.x,origin)->df.lolli

names(df.lolli) <- c(
  "Hugo_Symbol",
  "Protein_Change",
  "AA_Position",
  "Sample_ID",
  "Mutation_Type",
  "Chromosome",
  "Start_Position",
  "End_Position",
  "Reference_Allele",
  "Variant_Allele",
  "Center"
)
#write.table(df.lolli, file='test.tsv', quote=FALSE, sep='\t', col.names = TRUE,row.names=FALSE)

##plot with g3viz. Color coding by type of origin (WB/cfDNA), usually by Mutation_Type(SNV/stopgain/frameshift)
plot.options <- g3Lollipop.theme(theme.name = "cbioportal",
                                 title.text = "TP53",
                                 y.axis.label = "# of Mutations")

g3Lollipop(df.lolli%>%filter(Hugo_Symbol=="TP53"),
           gene.symbol = "TP53",
           btn.style = "gray", # gray-style chart download buttons
           plot.options = plot.options,
           factor.col = "Center",
           save.png.btn	= FALSE,
           save.svg.btn = FALSE,
           output.filename = "cbioportal_theme")

#####  Lolliplot for BRCA1/2 muts-----------------------------------------------------------------------
df.brca_germline%>%
  separate(.,AAChange,
           into=c("transcript1","rest"),
           sep=",",
           remove=TRUE,
           convert=FALSE)%>%
  separate(.,transcript1,
           into=c("gene","transcript_name","exon","DNAchange","amino_acid_change"),
           sep = ":",
           remove = TRUE,
           convert = FALSE)%>%
  #mutate(ExonicFunc = ifelse(ExonicFunc=="nonsynonymous SNV","Missense","Truncating"))%>%
  mutate(AA_pos = as.character(amino_acid_change))%>%
  separate(.,AA_pos,
           into=c("p","AA_pos1"),
           sep = "p.",
           remove = TRUE,
           convert = FALSE)%>%
  separate(.,AA_pos1,
           into=c("AA_pos2","rest"),
           sep = "fs",
           remove = TRUE,
           convert = FALSE)%>%
  mutate(AA_position = extract_numeric(AA_pos2))%>%
  dplyr::select(Gene,amino_acid_change,AA_position,Sample,ExonicFunc,Chr,Start,End,Ref,Alt)->df.lolli

names(df.lolli) <- c(
  "Hugo_Symbol",
  "Protein_Change",
  "AA_Position",
  "Sample_ID",
  "Mutation_Type",
  "Chromosome",
  "Start_Position",
  "End_Position",
  "Reference_Allele",
  "Variant_Allele"
)

##plot with g3viz
plot.options <- g3Lollipop.theme(theme.name = "cbioportal",
                                 title.text = "BRCA1",
                                 y.axis.label = "# of Mutations")

g3Lollipop(df.lolli%>%filter(Hugo_Symbol=="BRCA1"),
           gene.symbol = "BRCA1",
           btn.style = "gray", # gray-style chart download buttons
           plot.options = plot.options,
           factor.col = "Mutation_Type",
           save.png.btn	= FALSE,
           save.svg.btn = FALSE,
           output.filename = "cbioportal_theme")

##plot with g3viz
plot.options <- g3Lollipop.theme(theme.name = "cbioportal",
                                 title.text = "BRCA2",
                                 y.axis.label = "# of Mutations")

g3Lollipop(df.lolli%>%filter(Hugo_Symbol=="BRCA2"),
           gene.symbol = "BRCA2",
           btn.style = "gray", # gray-style chart download buttons
           plot.options = plot.options,
           factor.col = "Mutation_Type",
           save.png.btn	= FALSE,
           save.svg.btn = FALSE,
           output.filename = "cbioportal_theme")


#Confoundation by CH in HRD diagnostic ->ATM and CHEK2
ATMandCHEK2%>%
  filter(Gene == "ATM"| Gene == "CHEK2")%>%
  filter(tag == "true" | tag == "tumour")%>%
  mutate(origin = ifelse(tag == "true", "CH", 
                         ifelse(tag == "tumour", "tumour", NA)))->ATMandCHEK2

# Create a table of tag counts
tag_counts <- table(ATMandCHEK2$origin)

# Calculate the percentage of variants in each gene category
tag_percents <- prop.table(tag_counts) * 100

# Convert the table object to a data frame
tag_percents_df <- as.data.frame(tag_percents)

# Create a bar plot of the percentage of variants in each gene category
ggplot(tag_percents_df, aes(x = Var1, y = Freq, fill = )) +
  geom_bar(stat = "identity", fill = "#486081") +
  labs(title = "Percentage of Tumour or CH derived ATM and CHEK2 variants",
       x = "Origin",
       y = "Percentage")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.35, face = "italic", size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.35, face = "italic", size = 18),
        plot.title = element_text(size = 16, face = "bold"))->p.origin.hrd

png("output/figures/p.origin.hrd.png",width=10, height=10,units="in",res=500,type="cairo")
p.origin.hrd
dev.off()

#how many patients in total?
select(ATMandCHEK2,Patient.ID, origin)->hrd.diagnostic.pat
# calculate percentage and count of patients in each category
df_summary <- hrd.diagnostic.pat %>% 
  group_by(origin) %>% 
  summarize(count = n()) %>% 
  mutate(percentage = count / sum(count))

# create doughnut plot
ggplot(df_summary, aes(x = "", y = percentage, fill = origin)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(label = paste0(count, "\n(", round(percentage * 100), "%)")),
            position = position_stack(vjust = 0.5), 
            size = 10, color = "white") +
  scale_fill_manual(values = c("CH" = "#486081", "tumour" = "lightsteelblue")) +
  labs(title = "Total of 17 patients ",fill = "Origin")+
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.35, face = "italic", size = 16),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))->p.origin.hrd.total.pat
    
png("output/figures/p.origin.hrd.total.pat.png",width=10, height=10,units="in",res=500,type="cairo")
p.origin.hrd.total.pat
dev.off()

ATMandCHEK2_Cosmic <- ATMandCHEK2 %>% 
  mutate(COSMIC = case_when(
    str_detect(cosmic92_coding, "ovary|breast") ~ "Ovary",
    str_detect(cosmic92_coding, "lymphoid") ~ "CH",
    TRUE ~ NA_character_
  ))




#####  how many counts per sample ####
df.filtered%>%
  filter(Material == "cf")-> test
  
# Count the number of mutations for each sample
sample_counts <- test %>% group_by(Sample) %>% summarise(count = n())

# Count the number of samples with each count of mutations
no_counts <- sample_counts %>% group_by(count) %>% summarise(num_samples = n())

ggplot(no_counts, aes(x=count, y=num_samples))+
  geom_bar(stat = "identity", fill = "#486081") +
  geom_text(aes(label = count), vjust = -0.5)

sample_counts%>%
  filter(count >= 30)%>%
  mutate(Sample_orig = Sample)%>%
  left_join(.,ids,by = "Sample_orig")->test2 #candidates for repetition