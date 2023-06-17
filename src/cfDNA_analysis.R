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

##samples failed in sequencing
failedSamples <-c('OvCA_44_C1D1_cf','OvCA_45_C1D1_cf','OvCA_46_C1D1_cf','OvCA_48_C1D1_cf','OvCA_50_C1D1_cf','OvCA_54_C1D1_cf','OvCA_93_C1D1_cf',
                  'OvCA_11_C1D1_cf','OvCA_40_C1D1_cf','OvCA_53_C1D1_cf','OvCA_65_C1D1_cf')

##relevant variables
variables <- c("Patient.ID","Sample_orig","mutID","position","Sample", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "GeneDetail", "ExonicFunc", "AAChange", "cytoBand","readDepth", "TR1", "TR1_plus", "TR1_minus", "TR2", "TR2_plus", "TR2_minus", "TVAF", "AF", "avsnp150","cosmic92_coding","snp","mutFreq","p.binom","n.mut","n.material","sum_cf","sum_wb","Material","tag", "Patmut")

##interesting gene groups
ch_genes <- c("DNMT3A","TET2","ASXL1","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5")
tp53_genes <- c("TP53")
ppm1d_genes <- c("PPM1D")
brca_genes <- c("BRCA1","BRCA2")
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")
ch_genes_without_HRD <- c("DNMT3A","TET2","ASXL1","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5")
Categories<-c('CH','HRD','other','TP53')

#dataframe with mutation calls from cfDNA             
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

#####  Identity check via SNP  ####
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
            ylab = "whole-blood VAF", 
            xlab = "cfDNA VAF")->p.cfDNACor
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
  mutate(gene = ifelse(is.element(Gene.x,ch_genes_without_HRD),"CH",
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

#####  Plot that shows VAF WB vs VAF ctDNA including color for group of mutation - Filter 1 ####
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
  filter(p.binom.x <= -10) %>%
  filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  filter(ExonicFunc.x != "synonymous SNV")%>%
  filter(AF.x<0.1)%>%
  filter(snp.x==FALSE)%>%
  filter(TVAF.x>0.01|TVAF.y>0.01)%>%
  filter(TR2.y > 19|TR2.x>19)%>%
  ggplot(aes(x=TVAF.x,y=TVAF.y,
             color=tag.x,
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

png("output/figures/p.cf.corr.filter1.png",width=10, height=6,units="in",res=500,type="cairo")
p.cf.corr
dev.off()

##### Correlation plot Filter 1 - TP53 mutations only 
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(!is.element(Sample.x,failedSamples))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes_without_HRD),"CH",
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
  filter(Gene.x=="TP53")%>%
  #filter(cosmic_ovary)%>%
  ggplot(aes(x=TVAF.x,y=TVAF.y,color=cosmic_ovary))+
  geom_point(size=2)+
  geom_abline(slope=1)+
  scale_color_viridis(discrete=TRUE)+
  scale_x_log10(limits=c(0.0005,0.5)) +
  scale_y_log10(limits=c(0.0005,0.5)) +
  theme_minimal() -> p.TP53_cosmic_cf_wb
p.TP53_cosmic_cf_wb

png("output/figures/p.TP53_cosmic_cf_wb.png",width=10, height=6,units="in",res=500,type="cairo")
p.TP53_cosmic_cf_wb
dev.off()


#####  Max: Mutationspectrum for mutations in WB and cf only mutations by group------------------
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
  filter(tag.x != "false"|tag.x != "germline")%>%
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  dplyr::select(gene, compartment) %>% 
  table %>% 
  data.frame %>%
  ggplot(aes(x=reorder(gene, Freq), y=Freq, fill=compartment)) +
  geom_bar(stat="identity", width=0.6, position = "stack")+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(position = "right")+
  ylab("Gene Mutation Frequency") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35),
        axis.ticks.y = element_blank())+
  coord_flip() -> p.mutprev
p.mutprev

### Mutationspectrum cf vs wb by gene 
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(!is.element(Sample.x,failedSamples))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes_without_HRD),"CH",
                       ifelse(is.element(Gene.x,tp53_genes),"TP53",
                              ifelse(is.element(Gene.x,hrd_genes),"HRD",
                                     ifelse(is.element(Gene.x,brca_genes),"BRCA",
                                            ifelse(is.element(Gene.x,ppm1d_genes),"PPM1D","other"))))))%>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding.x,"ovary")) %>%
  filter(p.binom.x <= -Inf) %>%
  filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  filter(ExonicFunc.x != "synonymous SNV")%>%
  filter(AF.x<0.1)%>%
  filter(snp.x==FALSE)%>%
  filter(TVAF.x>0.01|TVAF.y>0.01)%>%
  filter(TR2.y > 19|TR2.x>19)%>%
  filter(TVAF.x<0.35&TVAF.y<0.35)%>% #no germline variants
  filter(tag.x != "false"|tag.x != "germline")%>%
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  dplyr::select(Gene.x, compartment) %>% 
  table %>% 
  data.frame %>%
  ggplot(aes(x=reorder(Gene.x, Freq), y=Freq, fill=compartment)) +
  geom_bar(stat="identity", width=0.6, position = "stack")+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(position = "right")+
  ylab("Gene Mutation Frequency") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35),
        axis.ticks.y = element_blank())+
  coord_flip() -> p.mutprev
p.mutprev


##Lolliplots for TP53 mutations

full_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(!is.element(Sample.x,failedSamples))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes_without_HRD),"CH",
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
  filter(Gene.x=="TP53")%>%
  mutate(Gene=Gene.x,
         AAChange = AAChange.x,
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

###Klara:----------df.filtered.cf soll im Endeffekt bei cf_variant filtering rauskommen, ich überarbeite das gerade noch, ich lösche das hier wieder, wenn es dann funktioniert
### below line composition analysis
df.filtered_cf%>%
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
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.35, size = 16),
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

#####  Klara: Gene Mutation Prevalence Plot (plots number of gene-x-mutated patients)  #####
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

#####  Klara: Determine Overlap CH and HRD in plasma samples ####
df.filtered%>%
  filter(Material=="cf") %>% 
  filter(Visite == "C1D1")%>%
  filter(is.na(replicate))%>%
  filter(TVAF >= 0.01 & TVAF <= 0.3)->df.filtered.cf

nop<-df.filtered.cf%>%
  select(.,Patient.ID)%>%
  unique()%>%nrow()
source("src/global_functions_themes.R")

df.filtered.cf %>% 
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


#####Klara:  Confounding by CH in HRD diagnostic ->ATM and CHEK2  ####
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
#####  Barplot cf components WBC-Tumour derived #####
load ('data/interim/workingdata.RDATA')

#canonical CH Gene from this paper: Clonal hematopoiesis detection in patients with cancer using cell-free DNA sequencing
canonicalCHgene<-c("TET2", "DNMT3A", "ASXL1", "JAK2", "PPM1D", "SF3B1")
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")
myeloidDKMS <-c("SRSF2","TP53", "U2AF1", "CBL", "IDH1", "IDH2", "BCOR", "BCORL1", "EZH2", "STAG2",
                "GNAS","GNB1","KRAS", "NRAS", "WT1", "MYD88", "STAT3", "CALR", "CEBPA", 
                "CSF3R", "ETV6", "FLT3","GATA2","GATA1","KIT","MPL","NPM1","PTPN11","RUNX1","SETBP1","NF1","PHF6")

workingdata%>%
  filter(visit_material == "C1D1_cf")%>%
  filter(tag == "true"| tag=="cf-only")%>%
  mutate(category = ifelse(is.element(Gene, canonicalCHgene), "canonical CH genes",
                           ifelse(is.element(Gene, myeloidDKMS), "other myeloid driver genes",
                                  ifelse(is.element(Gene, hrd_genes), "HRD genes", "other")))) ->cf.data
cf.data%>%
  nrow()->totalmutations

## Tumor or WB derived plot
cf.data %>%
  filter(tag == "true")%>%
  group_by(category) %>%
  summarize(count = n())%>%
  mutate(Percentage = count / totalmutations)%>%
  mutate(tag = paste("WBC derived"))-> plotdata1
  
cf.data %>%
  filter(tag == "cf-only")%>%
  group_by(category) %>%
  summarize(count = n())%>%
  mutate(Percentage = count / totalmutations)%>%
  mutate(tag = paste("Tumour derived"))-> plotdata2

bind_rows(plotdata1, plotdata2)%>%
  ggplot(aes(x = tag, y = Percentage, fill = category)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = " ", y = "Percentage of total mutations", fill = "Gene category") +
  scale_fill_manual(values=c("#deebf7", "#7A96B8", "#56779F", "#394F6A")) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))->p.percentageTumorWBC

png("output/figures/p.percentageTumorWBC.png",width=8, height=6,units="in",res=500,type="cairo")
p.percentageTumorWBC
dev.off()

## Tumor or WB derived with genes plot
#Prepare dataset
mutation_counts <- cf.data %>%
  group_by(Gene, tag) %>%
  tally()
mutation_counts_wide <- mutation_counts %>%
  spread(tag, n, fill = 0)

mutation_counts_wide <- mutation_counts_wide %>%
  mutate(Total = true + `cf-only`)

mutation_counts_summary <- mutation_counts_wide %>%
  pivot_longer(cols = c("true", "cf-only"), names_to = "origin", values_to = "Count") %>%
  mutate(Tag = ifelse(origin == "true", "true", "cf-only")) %>%
  group_by(Gene, origin) %>%
  summarize(Total = sum(Count)) %>%
  arrange(desc(Total))

# Create  barplot
ggplot(mutation_counts_summary, aes(x = reorder(Gene, -Total), y = Total, fill = origin)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#7A96B8", "#394F6A"), labels = c("WB derived", "Tumour derived")) +
  xlab("Gene") +
  ylab("Number of Mutations") +
  ggtitle("Mutation Counts by Gene according to origin") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))+
  scale_y_continuous(breaks = seq(0,65, 5))->p.GenecountTumorWBC

png("output/figures/p.GenecountTumorWBC.png",width=8, height=6,units="in",res=500,type="cairo")
p.GenecountTumorWBC
dev.off()

