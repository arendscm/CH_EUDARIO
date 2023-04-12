# ______________________________________________________________________________
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: from variantcalls to excel list of filtered results
#
# Input: seqdata
#
# Output: Excel list of filtered results
#
# ______________________________________________________________________________
########   Dependencies   #####
library(base)
library(dplyr)
library(xlsx)
library(stringr)
library(reshape)
library(tidyr)
library(readxl)
library(reshape2)


########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')

##load hotspot data
load('data/interim/hotspots.RData')

######## Get Patient ids
source("src/ids.R")
source("src/global_functions_themes.R")

##HRD genes
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3", "BRCA1", "BRCA2")
ch_genes <- c("DNMT3A", "TET2" ,  "JAK2" ,  "ASXL1" , "SF3B1" , "SRSF2" , "TP53"  , "U2AF1" , "PPM1D" , "CBL"  ,  "IDH1"  , "IDH2"  , "BCOR"  , "BCORL1", "EZH2" ,  "RAD21" , "STAG2" , "CHEK2" , "GNAS"  , "GNB1"  , "ATM"   , "KRAS" ,  "NRAS",   "WT1" ,   "MYD88" ,
              "STAT3" , "BRCC3" , "CALR"  , "CEBPA" , "CSF3R" , "ETV6"  , "FLT3" ,  "GATA2" , "GATA1" , "KIT" ,   "MPL" ,   "NPM1" ,  "PTPN11" ,"RUNX1" , "SETBP1" ,"NF1"  ,  "PHF6")


##identify variants with VAF < 1% through paired WB - plasma samples
##dataframe with c1d1 plasma calls
df %>% 
  filter(Material=="cf") %>% 
  filter(Visite == "C1D1")%>%
  filter(is.na(replicate))%>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))-> df.cf

#data frame with mutation calls from WB samples that have matched cfDNA samples
df %>% 
  filter(is.element(Patient.ID,df.cf$Patient.ID)) %>% 
  filter(is.na(replicate))%>%
  filter(Material=="wb") %>% 
  filter(Visite == "C1D1")%>%
  mutate(cfID=paste(Patient.ID,position,sep="_")) -> df.cf_wb

##identify rare clones that were called in wb and plasma
full_join(df.cf,df.cf_wb,by="cfID") %>% 
  #filter(!is.element(Sample.x,mismatch))%>%
  mutate(TVAF.y = ifelse(is.na(TVAF.y),0,TVAF.y)) %>% 
  mutate(TVAF.x = ifelse(is.na(TVAF.x),0,TVAF.x)) %>%
  filter(p.binom.y == -Inf) %>%
  filter(mutFreq.y < 0.2*n.lane.y) %>%
  filter(Func.x == "exonic"|Func.x == "splicing"|Func.x == "exonic;splicing") %>%
  filter(ExonicFunc.x != "synonymous SNV")%>%
  filter(AF.x<0.001)%>%
  filter(TVAF.x>0.001&TVAF.y>0.001)%>%
  filter(TR2.y > 9&TR2.x>9)%>%
  filter(TVAF.y < 0.35)-> df.rare

df.cf_wb$Patient.ID %>% unique %>% length -> nop

df.rare %>% 
  dplyr::select(Sample.y, Gene.y) %>% 
  data.frame %>% 
  unique %>% 
  dplyr::select(Gene.y) %>% 
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
  scale_y_continuous(labels = percent,limits=c(0,0.5), position = "right")+
  ylab("Gene Mutation Prevalence [%]") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  coord_flip() + 
  scale_fill_manual(values = c("non HRD" = "#486081", "HRD" = "#88acd4")) -> p.mutprev
p.mutprev

png("output/figures/mutprev.png",width=8, height=6,units="in",res=500,type="cairo")
p.mutprev
dev.off()



#####mutation frequency with exonic funciton
df.rare %>%
  dplyr::select(Gene.y,ExonicFunc.y) %>% 
  mutate(ExonicFunc.y = replace(ExonicFunc.y,ExonicFunc.y == ".","splice mutation")) %>%
  data.frame %>% 
  group_by(Gene.y) %>% 
  summarise(Gene.freq = n())%>%
  as.data.frame %>% 
  full_join(df.rare,.,by = "Gene.y") -> df.genefreq

df.genefreq %>%
  dplyr::select(Gene.y,ExonicFunc.y) %>% 
  mutate(ExonicFunc.y = replace(ExonicFunc.y,ExonicFunc.y == ".","splice mutation")) %>%
  data.frame %>% 
  table %>% 
  data.frame %>%
  ggplot(aes(x=reorder(Gene.y,Freq), y=Freq,fill=ExonicFunc.y)) +
  geom_bar(stat="identity", width=0.6)+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(position = "right")+
  ylab("Mutation Frequency") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  scale_fill_viridis(discrete=TRUE) +
  labs(fill="")+
  coord_flip() -> p.mutfreq
p.mutfreq


png("output/figures/mutfreq.png",width=8, height=6,units="in",res=500,type="cairo")
p.mutfreq
dev.off()

rm(p.mutfreq)
rm(df.genefreq)