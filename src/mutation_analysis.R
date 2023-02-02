# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: Mutation analysis and plots from filtered mutational data
#
# Input: df.filtered from data/interim/seqdata_filtered.RData
#
# Output: Excel list of filtered results, plots, ...
#
# ==============================================================================
########   Dependencies   ####
library(base)
library(dplyr)
library(stringr)
library(reshape)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(viridis)
library(reshape)
library(ggpubr)
library(g3viz)
library(maftools)

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata_filtered.RData')

######## Get Patient ids
source("src/ids.R")

######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")

########   Mutational Analysis------------------------------------------------------------------------------

##PLOTs

#Gene Mutation Prevalence Plot (plots number of gene-x-mutated patients)
nop <- ids%>%
  filter(Patient.ID!=0)%>%
  select(.,Patient.ID)%>%
  unique()%>%nrow #number of patients

hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")

df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
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
  scale_y_continuous(labels = percent,limits=c(0,0.3), position = "right")+
  ylab("Gene Mutation Prevalence [%]") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  coord_flip() + 
  scale_fill_manual(values = c("non HRD" = "#486081", "HRD" = "#88acd4")) -> p.mutprev
p.mutprev

png("output/figures/mutprev.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev
dev.off()

rm(prev.table)
rm(p.mutprev)

#CH prevalence by BRCA status
source("src/brca_germline.R")

df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  dplyr::select(Patient.ID) %>% unique -> id.ch
id.brca_germline %>% 
  mutate(CH = ifelse(is.element(Patient.ID,id.ch$Patient.ID),1,0))%>%
  dplyr::select(brca_germline,CH) %>% 
  table -> brca_status

##prevalence plot by BRCA status
#prevalences in BRCA wildtype patients

n.brcamut <- id.brca_germline %>% 
  mutate(CH = ifelse(is.element(Patient.ID,id.ch$Patient.ID),1,0))%>%
  dplyr::select(brca_germline) %>% sum 

df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  left_join(.,id.brca_germline,by = "Patient.ID")%>%
  filter(brca_germline==0)%>%
  dplyr::select(Sample, Gene) %>% 
  data.frame %>% 
  unique %>% 
  dplyr::select(Gene) %>% 
  table %>% 
  data.frame %>% 
  mutate(prev = Freq/(nop-n.brcamut)) %>% 
  arrange(prev) %>%
  mutate(brca = 0) -> prev.table_brca0
names(prev.table_brca0)<- c( "Gene","Freq","prev","brca")

#prevalences in brca mutated patients

df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  left_join(.,id.brca_germline,by = "Patient.ID")%>%
  filter(brca_germline==1)%>%
  dplyr::select(Sample, Gene) %>% 
  data.frame %>% 
  unique %>% 
  dplyr::select(Gene) %>% 
  table %>% 
  data.frame %>% 
  mutate(prev = Freq/n.brcamut) %>% 
  arrange(prev) %>%
  mutate(brca=1) -> prev.table_brca1
names(prev.table_brca1)<- c( "Gene","Freq","prev","brca")

full_join(prev.table_brca0,prev.table_brca1) %>% dplyr::select(Gene) %>% unique -> gene
gene %>% left_join(prev.table_brca0) %>% mutate(prev = ifelse(is.na(prev),0,prev),brca=0)->brca0
gene %>% left_join(prev.table_brca1)%>% mutate(prev = ifelse(is.na(prev),0,prev),brca=1)->brca1
full_join(brca0,brca1)->prev.brca

prev.brca  %>%
  ggplot(aes(x=reorder(Gene, prev), y=prev, fill = factor(brca),group=brca)) +
  geom_bar(stat="identity", width=0.6,position=position_dodge())+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(labels = percent,limits=c(0,0.4), position = "right")+
  ylab("Gene Mutation Prevalence [%]") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  scale_fill_viridis(discrete=TRUE) +
  scale_fill_manual(values = c("0" = "#486081", "1" = "#88acd4")) +
  coord_flip() -> p.mutprev
p.mutprev

png("plots/mutprev-BRCA.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev
dev.off()

########   GenVisar - waterfall plot #### 

library(viridisLite)

library(GenVisR)

###turn df.filtered into MAF compatible format
df.maf <- makeMAF(df.filtered.c1d1%>% filter(tag=="true"))


# count sample number
n <- length(unique(df.maf$Tumor_Sample_Barcode))
# input x colors for x types of mutations
#color_R <- c("firebrick4", 
#"darkgoldenrod1",
#"orangered3", 
#"deepskyblue4",
#"darkgreen")
            
# modify Layers of GenVisR::waterfall to look pretty
# sampRecurLayer: ggplot Layer to be added to the left subplot
sampRecurLayer <- theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        panel.grid = element_blank(),
                        plot.margin = margin(t = 5, r = 5, b = 2, l = 5, unit = "pt"),
                        axis.text.x = element_text(color = "black"))

# mainLayer: ggplot layer to be added to the main waterfall plot
mainLayer <- theme(axis.text.y = element_text(size = 10.5, color = "black", face = "italic"),
                   axis.text.x = element_text(size = 12, color = "black", face = "italic"),
                   plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"),
                   legend.title = element_text(size = 12),
                   axis.title.x = element_blank())

gene.order <- genes$Hugo_Symbol                      ## vector specifying target genes and gene order

variants <- unique(OvCAmaf$Variant_Classification)    ## vector specifying types of mutation

# or
# variants <- c("Nonsynonymous SNV",                 ## you can choose any type of mutation you want
#              "Frameshift Deletion",                ## for the "costum" file type
#             "Splice Site", 
#              "Stop-Gain"
#               )

names(OvCAmaf) <- c("sample", "gene", "variant_class")  ## rename dataframe to fit costum type

# Generate plot and save as .pdf
#with MAF format
waterfall(df.maf, 
          fileType = "MAF",
          rmvSilent = FALSE,               ## don't show silent mutations (we don't have any)
          mainDropMut = TRUE,             ## unused mutation types will be dropped from the legend
          mainPalette = viridis(5),          ## previously specified colors
          plotMutBurden = FALSE,          ## mutation burden plot at the top
          mainLayer = mainLayer,          ## specify defined layers
          sampRecurLayer = sampRecurLayer, 
          #geneOrder = gene.order,         ## order to plot the genes (default: decreasing gene freq.)
          #variant_class_order = Mutation_Type, 
          section_heights = c(0.15, 1),   ## relative size of the different plots
          mainXlabel = FALSE,
          mainLabelSize = 0,
          main_geneLabSize = 1)


########   number of mutations plot------------------------------------------------------------
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%  
  filter(TVAF >= 0.01) %>%
  dplyr::select(Gene,ExonicFunc) %>% 
  mutate(ExonicFunc = replace(ExonicFunc,ExonicFunc == ".","splice mutation")) %>%
  data.frame %>% 
  group_by(Gene) %>% 
  summarise(Gene.freq = n())%>%
  as.data.frame %>% 
  full_join(df.filtered %>% filter(tag == "true"),.,by = "Gene") -> df.genefreq

##mutation frequency according to type of mutation 
df.genefreq %>%
  dplyr::select(Gene,ExonicFunc) %>% 
  mutate(ExonicFunc = replace(ExonicFunc,ExonicFunc == ".","splice mutation")) %>%
  data.frame %>% 
  table %>% 
  data.frame %>%
  ggplot(aes(x=reorder(Gene,Freq), y=Freq,fill=ExonicFunc)) +
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


png("output/figures/mutfreq.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutfreq
dev.off()

rm(p.mutfreq)
rm(df.genefreq)



########   SERIAL SAMPLES ####

#find Clones found in C1D1 and EOT for each patient
df.filtered%>%
  filter(Visite== "EOT")%>%
  select(.,Patient.ID)%>%
  filter(Patient.ID != "0")%>%
  unique()%>%
  unlist()->EOT_ids


#für ersten Durchlauf -> table zur Verfügung stellen
Dynamicallpat<-Dynamic

for(current_patient in EOT_ids)
{print(current_patient)
  #Datensätze auf Pat anpassen 
  df%>%
    filter(Patient.ID==current_patient)%>%
    filter(mutFreq<=15)%>%
    filter(snp=="FALSE")%>%
    filter(Visite != "cf-C1D1")->df.pat
  df.pat%>%
    filter(Visite=="C1D1")->df.c1d1pat
  df.pat%>%
    filter(Visite=="EOT")->df.eotpat
  
  #Überschneidungen finden
  inner_join(df.c1d1pat,df.eotpat,by="Patmut")->C1D1EOTdynamic
  
  #Matching Mutationen markieren 
  C1D1EOTdynamic%>%
    mutate(match="yes")->C1D1EOTdynamic
  select(C1D1EOTdynamic,Patmut,match)->Dynamicmut
  
  #mutationen aus dem großen pool patpool rausholen 
  left_join(df.pat,Dynamicmut,by="Patmut")->Dynamic
  Dynamic%>%
    filter(match=="yes")->Dynamic
  
  ##bis jetzt sind einfach alle matching mutationen mit angewandten Filtern drin, viele sind aber sehr klein!
  ##jetzt: C1D1 >0,5 filtern und dann wieder im EOT pool gucken, ob es matches gibt!
  Dynamic%>%
    filter(Visite=="C1D1")%>%
    filter(TVAF >= 0.005)->DynamicC
  Dynamic%>%
    filter(Visite=="EOT")->DynamicE
  
  #Überschneidungen finden
  inner_join(DynamicC,DynamicE,by="Patmut")->C1D1EOTdynamic
  
  #Matching Mutationen markieren 
  C1D1EOTdynamic%>%
    mutate(match2="yes")->C1D1EOTdynamic
  select(C1D1EOTdynamic,Patmut,match2)->Dynamicmut
  
  #mutationen aus dem großen pool patpool rausholen 
  left_join(Dynamic,Dynamicmut,by="Patmut")->Dynamicx
  Dynamicx%>%
    filter(match2=="yes")->Dynamicx
  
  ##Jetzt müssen noch EOT >0.5 gerettet werden die in C1D1 noch kein match hatten!
  DynamicE%>%
    filter(TVAF >=0.005)->DynamicEresc
  #die, die noch kein match haben
  anti_join(DynamicEresc,DynamicC, by="Patmut")->DynamicEresc
  #
  DynamicEresc%>%
    mutate(match3="yes")->DynamicEresc
  select(DynamicEresc,Patmut,match3)->Dynamicmut
  # aus C1D1 Pool rescuen
  semi_join(df.c1d1pat,Dynamicmut, by="Patmut")->DynamicCresc
  bind_rows(DynamicEresc,DynamicCresc)->Dynamicy
  
  #Matches und rescued zusammenführen
  bind_rows(Dynamicx,Dynamicy)->Dynamic
  
  #Bind Dynamic with other Pat Dynamics files
  bind_rows(Dynamicallpat,Dynamic)%>%
    unique()%>%
    filter(cf != "1")->Dynamicallpat
}


rm(C1D1EOTdynamic)+
rm(df.pat)+
rm(DynamicC)+
rm(DynamicCresc)+
rm(DynamicE)+
rm(DynamicEresc)+
rm(Dynamicmut)+
rm(Dynamicx)+
rm(Dynamicy)+
rm(Dynamic)



#write to excel file
filename="C1D1EOT/Dynamicallpatients.xlsx"
#write.xlsx(Dynamicallpat,filename, sheetName = "Dynamic",append=TRUE)

#Create Dynamic plots for each EOT ID
#setwd to where images should be saved to!
setwd("H:/Meine Ablage/C1D1EOT/Dynamicplot")
for (current_patient in EOT_ids)
{print(current_patient)
  Dynamicallpat%>%
    filter(Patient.ID == current_patient)->Dynamic
  ###Plot generieren
  Dynamic%>%
    ggplot(aes(y=TVAF, x=Visite, colour=Gene, group=Patmut))+
    geom_point(size=5,alpha=0.3)+
    geom_line(size=1)+
    theme_minimal()+
    scale_y_continuous(limits=c(0,0.125))+
    labs(title="Clone Dynamics")->p.C1D1EOTdynamicpat1
  ### create image file
  png(paste0(current_patient,"_C1D1EOT.png"),
      width=10,
      height=6,
      units="in",
      res=500,
      type="cairo")
  ### plot image to file
  print(p.C1D1EOTdynamicpat1)
  ### close file again
  dev.off()
}

rm(p.C1D1EOTdynamicpat1)+
rm(current_patient)+
rm(EOT_ids)+
rm(Dynamic)+
rm(Dynamicallpat)
rm(filename)

setwd("H:/Meine Ablage")
####  identity check via SNPs ####
df %>% 
  filter(Visite=="EOT") %>% 
  mutate(ID=paste(Patient.ID,position))-> df.eot1

ids %>% filter(Visite=="EOT") -> eotsamples

df %>% 
  filter(Visite=="C1D1") %>% 
  filter(is.element(Patient.ID,eotsamples$Patient.ID)) %>%
  mutate(ID=paste(Patient.ID,position)) -> df.c1d0

left_join(df.eot1,df.c1d0,by="ID") %>% 
  filter(snp.x == 1) %>% 
  ggplot(aes(x=factor(Patient.ID.x),y=TVAF.x-TVAF.y)) +
  geom_point()

##serial samples dynamics plot
df %>% filter(is.element(Patient.ID,eotsamples$Patient.ID)) %>% filter(cf==0)-> df.eot

df.eot %>% 
  filter(serial.mut>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(mutFreq < 10) %>% 
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.008) %>%
  filter(TVAF < 0.38) %>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=position,color=Gene),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ Patient.ID, ncol=6, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.serial

png("p.serial.png",width=6, height=4,units="in",res=500,type="cairo")
p.serial
dev.off()

###Serial samples by brca status (question: do dynamics unter PARP Inhb. differ depending on BRCA status?)
df.eot %>% 
  left_join(.,id.brca_germline,by = "Patient.ID")%>%
  filter(Gene=="PPM1D"|Gene=="TP53")%>%
  filter(serial.mut>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(mutFreq < 10) %>% 
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.008) %>%
  filter(TVAF < 0.38) %>%
  ggplot() + 
  geom_point(aes(x=Visite,y=TVAF,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=Visite,y=TVAF,group=position,color=Gene),size=1*1,na.rm=FALSE) + 
  facet_wrap(~ brca_germline, ncol=2, scales="free", dir="h") +
  scale_y_continuous(limits = c(0,0.26)) +
  labs(x="Time in days",y="Variant allele frequency",colour="Mutated Gene") +
  theme_minimal()-> p.serial

##TEST SERIAL SAMPLES relative >- for this we need the timedifference between d1 and eot

df.eot %>% 
  filter(serial.mut>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(mutFreq < 10) %>% 
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.005) %>%
  filter(TVAF < 0.37) %>%
  filter(Visite == "C1D1") %>%
  mutate(vaf_d1=TVAF) %>%
dplyr::select(Patient.ID,Gene,AAChange,position,vaf_d1)-> df.eotd1

df.eot %>% 
  filter(serial.mut>1)%>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  filter(AF<0.1)%>%
  filter(snp==FALSE)%>%
  filter(mutFreq < 10) %>% 
  group_by(Patient.ID,position) %>%
  mutate(maxVAF = max(TVAF)) %>%
  data.frame()%>%
  filter(maxVAF > 0.005) %>%
  filter(TVAF < 0.37) %>%
  filter(Visite == "EOT") %>%
  mutate(vaf_eot=TVAF) %>%
  dplyr::select(Patient.ID,Gene,AAChange,position,vaf_eot)-> df.eoteot

df.eot_rel <- full_join(df.eotd1,df.eoteot) %>% mutate(relvaf1 = vaf_d1/vaf_d1,
                                                       relvaf2 = vaf_eot/vaf_d1) %>%
  melt.data.frame(measure.vars = c("relvaf1","relvaf2"))

df.eot_rel %>%  
  ggplot() + 
  geom_point(aes(x=variable,y=value,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  geom_line(aes(x=variable,y=value,group=position,color=Gene),size=1*1,na.rm=FALSE) + 
  scale_y_log10() +
  labs(x="Timepoint",y="log(VAF change)",colour="Mutated Gene") +
  theme_graphicalabstract()-> p.serial

df.eot_rel %>% 
  filter(variable == "relvaf2") %>% 
  ggplot() + 
  geom_point(aes(x=Gene,y=value,color=Gene,group=Patient.ID),size=1.5,na.rm=FALSE) + 
  scale_y_log10() +
  labs(x="Gene",y="log(VAF change)",colour="Mutated Gene") +
  theme_graphicalabstract()-> p.serial

df.eot_rel %>% 
  filter(variable == "relvaf2") %>% 
  filter(is.element(Gene,c("CHEK2","PPM1D","DNMT3A","TP53","TET2","ASXL1")))%>%
  ggboxplot(., 
                 x = "Gene",
                 y = "value",
                 #facet.by = "variable",
                 # panel.labs = list(CHIP = c("positive","negative"), variable=c("Troponin","VCAM","hsCRP")), 
                 combine = TRUE,
                 color = "Gene", 
                 #palette = viridis(5),
                 xlab = "Gene",
                 ylab = "Growthrate",
                 title = "",
                 width = 0.3,
                 ylim = c(0,15),
                 size=0.8,
                 alpha=1,
                 repel=TRUE,
                 #yscale = "log10",
                 scales = "free",
                 add = c("jitter"))+
  theme_graphicalabstract() + 
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none",
        axis.title.y = element_text(face ="plain"),
        plot.title = element_text(hjust=0,face ="plain")) ->p.serial


########   cf DNA analysis ####
#nur Mutationen in cf
df%>%
  filter(Visite == "cf-C1D1")->df.cf
df%>%
  filter(Visite=="C1D1")->df.c1

anti_join(df.cf,df.c1,by="Patmut")->test
test%>%
  filter(Sample!="cf2-A7")%>%
  filter(TVAF>0.01 & TVAF <0.35)%>%
  filter(snp=="FALSE")%>%
  filter(mutFreq<10)%>%
  filter(FisherScore>1 & FisherScore <20)%>%
  filter(Func == "exonic"|Func=="splicing")->cfcalls

filename="cfcalls.xlsx"
write.xlsx(cfcalls,filename,append=TRUE)


#cf Mutations deren VAF x mal höher als in WB sind
df %>% filter(Visite == "cf-C1D1") %>% dplyr::select(ID) %>% unique -> cfsamples

#nur welche die auch cf haben und 4202008 raus da cf-2-A7 hunderte Mutationen hat!
df %>%
  filter(is.element(ID,cfsamples$ID)) %>%
  filter(ID != "4202008")%>%
  filter(TVAF <= 0.35 & TVAF >= 0.005)%>%
  filter(mutFreq <=10)->df.cf

select(df.cf,Chr,Start,End,Ref,Alt,Gene,Func,Visite,TVAF,AF,avsnp150,readDepth,FisherScore,EBScore,run, ID,Patmut,mutFreq)->df.cf

df.cf%>%
  filter(Visite=="C1D1")->WB
df.cf%>%
  filter(Visite=="cf-C1D1")->cf
inner_join(WB,cf,by="Patmut")%>%
  mutate(factor =(TVAF.y/TVAF.x))->WBcf
#.x=WB
#.y=cf

#filtering
WBcf%>%
  filter(TVAF.y >= (1.5*TVAF.x))->WBcffiltered 
#oder so, geht beides
WBcf%>%
  filter(TVAF.y/TVAF.x >1.5) ->WBcffiltered 


#Correlation Plot
WBcf%>% 
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

##interesting gene groups
variables <- c("Patient.ID","Sample_orig","cf","mutID","position","Sample", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "GeneDetail", "ExonicFunc", "AAChange", "cytoBand","readDepth", "TR1", "TR1_plus", "TR1_minus", "TR2", "TR2_plus", "TR2_minus", "TVAF", "AF", "avsnp150", "cosmic92_coding","snp","mutFreq")
ch_genes <- c("DNMT3A","TET2","ASXL1","PPM1D","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5")
tp53_genes <- c("TP53")
brca_genes <- c("BRCA1","BRCA2")
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")

#data frame with mutation calls from cfDNA             
df %>% 
  filter(cf==1) %>% 
  dplyr::select(variables,p.binom) %>%
  mutate(cfID=paste(Patient.ID,position))-> df.cf
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

##Rose Chart
# Create dataset
df.filtered.c1d1%>%
  filter(tag== "true")%>%
  filter(TVAF >= 0.01)%>%
  select(.,mutID,Gene,TVAF)%>%
  filter(is.element(Gene,c("CHEK2","PPM1D","DNMT3A","TP53","TET2","ATM")))->data
#data = data %>% arrange(Gene, TVAF)
#orders mutID acording to TVAF in each Gene group

genes <-c("PPM1D","DNMT3A","CHEK2","TP53","TET2","ATM")

# Initialize an empty list to store the filtered data for each gene
data.space<- list()

# Loop through each gene
for (gene in genes) {
  # Filter the data for the current gene
  filtered_data <- data[data$Gene == gene, ]
  
  # Add 4 empty rows to the bottom of the filtered data
  filtered_data <- rbind(filtered_data, data.frame(mutID = rep("", 4), Gene = rep("", 4), TVAF = rep(NA, 4)))
  
  # Add the filtered data to the list
  data.space[[gene]] <- filtered_data
}

# Join all the filtered data sets into one big data frame
final_data <- do.call(rbind, data.space)

#add id countin 1-...
final_data$id <- seq(1, nrow(final_data))

final_data->data

empty_bar <- 4


# prepare a data frame for base lines
base_data <- data %>% 
  group_by(Gene) %>% 
  summarize(start=min(id), end=max(id)) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))%>%
  filter(Gene != "")

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


# Make the plot
p.rosechart <- ggplot(final_data, aes(x=as.factor(id), y=TVAF, fill=Gene)) +      
  geom_bar(stat="identity", width=0.9) +
  
  
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.05, xend = start, yend = 0.05), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.1, xend = start, yend = 0.1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.2, xend = start, yend = 0.2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  annotate("text", x = rep(max(data$id),4), y = c(0,0.05,0.1,0.2), label = c("0","5%","10%","20%") , color="grey", size=1.5 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(stat="identity", width=0.9) +
  
  ylim(-0.1,0.35) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank())+
  coord_polar()+

  geom_segment(data=base_data, aes(x = start, y = -0.01, xend = end, yend = -0.01), colour = "grey", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = -0.032, label=Gene), colour = "grey", alpha=0.8, size=1.7, fontface="bold", inherit.aes = FALSE)


p.rosechart

png("output/figures/p.rosechart.png",width=8, height=8,units="in",res=500,type="cairo")
p.rosechart
dev.off()


