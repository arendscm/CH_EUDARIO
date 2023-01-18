# ==============================================================================
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
# ==============================================================================
########   Dependencies   #####
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


########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')

######## Get Patient ids
source("src/ids.R")

######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")

########   cf DNA analysis ####

##interesting gene groups
variables <- c("Patient.ID","Sample_orig","cf","mutID","position","Sample", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "GeneDetail", "ExonicFunc", "AAChange", "cytoBand","readDepth", "TR1", "TR1_plus", "TR1_minus", "TR2", "TR2_plus", "TR2_minus", "TVAF", "AF", "avsnp150", "cosmic92_coding","snp","mutFreq")
ch_genes <- c("DNMT3A","TET2","ASXL1","PPM1D","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5")
tp53_genes <- c("TP53")
brca_genes <- c("BRCA1","BRCA2")
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")

#data frame with mutation calls from cfDNA             
df %>% 
  filter(cf==1) %>% 
  filter(Visite == "C1D1")%>%
  dplyr::select(variables,p.binom) %>%
  mutate(cfID=paste(Patient.ID,position))-> df.cf
df %>% 
  filter(cf==0) %>% 
  filter(Visite == "C1D1")%>%
  dplyr::select(variables,p.binom) %>%
  mutate(cfID=paste(Patient.ID,position))-> df.wb

anti_join(df.cf, df.wb, by= "cfID")->df.cf_only
save(df.cf_only,file="data/interim/df.cf_only.RDATA")

#data frame with mutation calls from WB samples that have matched cfDNA samples
df %>% 
  filter(is.element(Sample,df.cf$Sample)) %>% 
  filter(cf==0) %>% 
  dplyr::select(variables,p.binom) %>%
  mutate(cfID=paste(Patient.ID,position)) -> df.cf_wb

##identity check via SNPs
left_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(snp.x == 1) %>% 
  ggplot(aes(x=Sample.x,y=TVAF.x-TVAF.y)) +
  geom_point()

mismatch <- c("2-E2","2-E8","2-H8","2-A7")

#Correlation Plot WB vs cfDNA
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(!is.element(Sample.x,mismatch))%>% 
  ggscatter(., 
            x = "TVAF.x", 
            y = "TVAF.y", 
            add = "reg.line", 
            #conf.int = TRUE, 
            cor.coef = TRUE, 
            cor.method = "pearson",
            size=1,
            xlab = "VAF whole blood", 
            ylab = "VAF plasma")->p.cfDNACor
p.cfDNACor


##Plot that shows VAF WB vs VAF ctDNA including color for group of mutation
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  filter(!is.element(Sample.x,mismatch))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  mutate(gene = ifelse(is.element(Gene.x,ch_genes),"CH",
                       ifelse(is.element(Gene.x,tp53_genes),"TP53",
                              ifelse(is.element(Gene.x,hrd_genes),"HRD",
                                     ifelse(is.element(Gene.x,brca_genes),"BRCA","other")))))%>%
  #filter(gene != "other") %>%
  mutate(cosmic_ovary = str_detect(cosmic92_coding.x,"ovary")) %>%
  filter(p.binom.x <= -10) %>%
  filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  filter(ExonicFunc.x != "synonymous SNV")%>%
  filter(AF.x<0.1)%>%
  filter(snp.x==FALSE)%>%
  filter(TVAF.x>0.005|TVAF.y>0.005)%>%
  filter(TR2.y > 14|TR2.x>14)%>%
  ggplot(aes(x=TVAF.x,y=TVAF.y,
             color=gene,
             #shape=ExonicFunc.x
             )) +
  geom_point(size=3)+
  geom_abline(slope=1,size=1,linetype=2,alpha=0.5)+
  scale_y_log10()+
  scale_x_log10()+
  #facet_wrap(~ Patient.ID.x, ncol=4, dir="h")+
  scale_color_viridis(discrete=TRUE)+
  ylab("whole blood VAF")+
  xlab("cfDNA VAF")+
  theme_minimal()

##TP53 mutations only 
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
  filter(p.binom.x == -Inf|p.binom.y == -Inf) %>%
  filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  filter(ExonicFunc.x != "synonymous SNV")%>%
  filter(AF.x<0.1)%>%
  #filter(snp.x==FALSE)%>%
  filter(TR2.x>9|TR2.y>9) %>% 
  filter(TVAF.x>0.005|TVAF.y>0.005)%>%
  filter(gene=="TP53")%>%
  ggplot(aes(x=TVAF.x,y=TVAF.y,color=gene))+
  geom_point(size=4)+
  geom_abline(slope=1)+
  scale_color_viridis(discrete=TRUE)+
  scale_x_log10(limits=c(0.0005,0.5)) +
  scale_y_log10(limits=c(0.0005,0.5)) +
  theme_minimal()

##detecting BRCA mutations in cfDNA (this will later on also be important when looking for BRCA reversion mutations)
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

########   BRCA and CH Status ####
##gehört hier eigentlich nicht hin zu cf Analysis...
#Ich denke wir sollten das in ein neuse Skript machen, was später durch TableOne ersetzt/erweitert wird
#Create table with BRCA and CH Status
for (id in ID)
{#####CH status
  print(id)
  df.filtered%>%
    filter(TVAF >= 0.005)%>%
    filter(ID == id)%>%
    filter(tag=="true")->CHposresults
  #nur CH Mutationen
  semi_join(CHposresults,CHgenes,by="Gene")->CHposresults
  nrow(CHposresults)->a
  if (a >= 1)
  {b<-"1"}
  if (a == 0)
  {b<-"0"}
  
  ####BRCA status
  df.BRCA%>%
    filter(ID== id)->BRCApos
  nrow(BRCApos)->c
  if (c >= 1)
  {d<-"1"}
  if (c == 0)
  {d<-"0"}
  
  tab2 <- matrix(c(id,b,a,d,c),ncol=5, byrow=TRUE)
  colnames(tab) <- c('Patient.ID','CH','nrow(CH)','BRCA','nrow(BRCA)')
  tab2 <- as.table(tab2)
  rbind(tab,tab2)->tab
  rm(a)
  rm(b)
  rm(c)
  rm(d)
  rm(id)
}

tab <- matrix(c(0,0,0,0,0),ncol=5, byrow=TRUE)
colnames(tab) <- c('ID','CH','nrow(CH)','BRCA','nrow(BRCA)')
tab <- as.table(tab)

tbl_df(tab)->BRCAStatus
BRCAStatus%>%
  filter(CH != 0)->BRCAStatus
#--> BRCAStatus = Table with BRCA and CH Status

#Bind BRCA Status with filtered results
as.numeric(BRCAStatus$ID)->BRCAStatus$ID
right_join(BRCAStatus,filtered_results_OvCA_fix,by="ID")->df.filtered

filteredmitBRCAStatus%>%
  filter(tag=="true")%>%
  filter(Visite=="C1D1")%>%
  filter(CH =="1")%>%
  filter(is.element(Gene,c("CHEK2","PPM1D","DNMT3A","TP53","TET2","ASXL1","ATM")))->a

#Create Boxplot according to BRCA Status
filteredmitBRCAStatus %>%
  ggplot(., aes(x=Gene,y=TVAF, fill=BRCA)) + 
  geom_boxplot()->p.cfboxplot1

filteredmitBRCAStatus%>%
  filter(is.element(Gene,c("PPM1D","DNMT3A","TP53","TET2")))->a

a %>%
  ggplot(., aes(x=Gene,y=TVAF, fill=BRCA)) + 
  geom_boxplot()->p.cfboxplot1

########   Lolliplot for TP53 muts------------------------------------------------------------------------
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

########   Lolliplot for BRCA1/2 muts-----------------------------------------------------------------------
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
