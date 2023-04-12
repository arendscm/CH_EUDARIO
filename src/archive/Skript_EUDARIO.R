# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: from variantcalls to excel list of filtered results and analysis
#
# Input: variantcalls as csv file
#
# Output: Excel list of filtered results, plots, ...
#
# ==============================================================================
########   Dependencies   #####
library(base)+
library(dplyr)+
library(ggplot2)+
library(xlsx)+
library(stringr)+
library(ggthemes)+
library(viridis)+
library(reshape)+
library(ggpubr)+
library(g3viz)+
library(tidyr)+
library(readxl)+
library(reshape2)

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########   QC sequencing runs  ####
HSM_1346 <- read_excel("QC/QC_runs.xlsx", sheet = "HSMetrics-P1346", 
                      col_types = c("text", "numeric", "numeric"))
GS_1346 <- read_excel("QC/QC_runs.xlsx", sheet = "General Statistics-P1346", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))
HSM_1519 <- read_excel("QC/QC_runs.xlsx", sheet = "HSMetrics-P1519", 
                       col_types = c("text", "numeric", "numeric"))
GS_1519 <- read_excel("QC/QC_runs.xlsx", sheet = "General Statistics-P1346", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))
full_join(HSM_1346,HSM_1519)->HSM
rm(HSM_1346)
rm(HSM_1519)
bind_rows(GS_1346,GS_1519)->GS
rm(GS_1346)
rm(GS_1519)

left_join(HSM,GS,by="Sample Name")->QC
rm(HSM)
rm(GS)
summary(QC)->d
rm(QC)

filename="QC/QC-sequencing runs_means.xlsx"
#write.xlsx(d,filename,append=TRUE)
rm(d)


########   Sample Registry and Int-Ext Pat ID ####

#Sample Registry -> what timepoint/sample do we have from each patient
Ovarial_Ca_LibPrep <- read_excel("data/external/Ovarial-Ca_LibPrep.xlsx", 
                                 sheet = "Pooling", skip=1)
select(Ovarial_Ca_LibPrep,'External Pat ID', 'Visite','Internal Pat ID', 'External Sample ID')->pool

pool%>%
  filter(Visite=="C1D1")->C1D1
pool%>%
  filter(Visite=="EOT")->EOT
pool%>%
  filter(Visite=="cf-C1D1")->cfC1D1
pool%>%
  filter(Visite=="cf-C7D1")->cfC7D1
pool%>%
  filter(Visite=="cf-EOT")->cfEOT

full_join(C1D1,EOT, by='External Pat ID')->list
full_join(list, cfC1D1, by= 'External Pat ID')->list
full_join(list, cfC7D1, by= 'External Pat ID')->list
full_join(list, cfEOT, by= 'External Pat ID')->list

filename="data/interim/Sample Registry.xlsx"
#write.xlsx(list, filename, sheetName = "list", append=TRUE)
rm(C1D1)
rm(cfC1D1)
rm(cfC7D1)
rm(cfEOT)
rm(EOT)
rm(pool)
rm(list)



#combine Ext and Int PatID
#komb aus Extraction Plan 
komb <- read_excel("Sample Registry/DNA Extraction Plan.xlsx", 
                   sheet = "Kombinieren")
#IntExt aus sample registry
IntExt <- read_excel("Sample Registry/Sample Registry.xlsx", 
                     sheet = "Int-Ext ID")
left_join(komb,IntExt, by='External Pat ID')->komb2

filename="Sample Registry/IntExtPatID.xlsx"
#write.xlsx(komb2, filename, sheetName = "list", append=TRUE)
rm(komb2)
rm(IntExt)
rm(komb)

########   Themes and Functions-----------------------------------------------------------------------------------
theme_graphicalabstract <- function(base_size=12, base_family="Roboto") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing.y = unit(0.2, "cm"),
           legend.title = element_text(),
           plot.margin=unit(c(10,5,5,5),"mm"),
           panel.spacing.x = unit(1,"line"),
           strip.background=element_rect(colour=NA,fill=NA)
           #strip.text = element_text(face="bold")
   ))
  
}

#Percent Function
percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "")
}

########   Data preparation:load variant calling table, Patient IDs, BRCA and tags-----------------------------------------------------------------
data1 <- read.table('variantcalls_P1346.csv',
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data2 <- read.table('variantcalls_P1519.csv',
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data1%>%
  mutate(run = "P1346")->data1
data2%>%
  mutate(run = "P1519")->data2 

#or load Variantcalls.RDATA from WD

##Patient ID table that identifies Sample IDs with Patient ID and timepoints
ids <- read_excel("Sample Registry/Ovarial-Ca_LibPrep.xlsx", 
                  sheet = "PatIDfix")
ids %>% filter(!is.na(Patient.ID)) -> ids 

# BRCA Exchange database for BRCA germline status
brcaexchange <- read.table("BRCA/BRCA_Exchange_Liste_shortend.csv",sep=";",header=TRUE)
#brcaexchange <- read.table("brca_exchange.tsv",sep="\t",header=TRUE)

# load tags 
tags <- read_excel("filtered_results.xlsx", 
                   sheet = "tags")

#internal Sample IDs
Sample_IDs <- read_excel("Sample Registry/Sample Registry.xlsx", 
                         sheet = "Sample IDs")
select(Sample_IDs,'Sample_orig','Sample.ID', 'External Sample ID', 'Visite', 'Internal Pat ID')->Sample_IDs
left_join(ids,Sample_IDs)->ids

rm(Sample_IDs)



########   Data Preparation: Restructure data for filtering---------------------------------------------------------------------
##relevant variables
variables <- c("Sample", "Patient.ID", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Func", "GeneDetail", "ExonicFunc", "AAChange", "cytoBand","readDepth", "TR1", "TR1_plus", "TR1_minus", "TR2", "TR2_plus", "TR2_minus", "TVAF", "AF", "avsnp150", "cosmic92_coding")

#hotspots
hotspots <- c("chr2_25234374_25234374",     #DNMT3A R882C
              "chr2_25234373_25234373",     #DNMT3A R882H
              "chr2_208248389_208248389",   #IDH1 R132
              "chr1_1815790_1815790",       #GNB1 K57
              "chr1_1815789_1815789",       #GNB1 K57
              "chr12_25245350_25245350",    #KRAS G12D
              "chr2_197402636_197402636",   #SF3B1 K666
              "chr2 197402635 197402635",   #SF3B1 K666
              "chr2_197402110_197402110",   #SF3B1 K700
              "chr17_76736877_76736877",    #SRSF2 P95
              "chr9_5073770_5073770",       #JAK2 V617
              "chr20_58909366_58909366"     #GNAS R844
)

##calculate mutation frequencies run-wise
for (file in c("data1", "data2")){
  print(file)
  eval(as.name(file)) %>%
    mutate(Sample_orig = Sample)%>%
    mutate(Sample = str_remove(Sample,"cf"))%>% #remove cf from Sample name and create column cf 1/0 instead
    mutate(cf = ifelse(str_detect(Sample_orig,"cf"),1,0))%>%
    mutate(AF = replace(AF, AF == ".", 0)) %>% 
    mutate(AF = as.numeric(AF)) %>%
    left_join(.,ids,by = "Sample_orig") %>% ##fuse with SampleID Table
    mutate(position = paste(Chr,Start,End,Ref,Alt,sep="_"))%>%  #mutation position with basechange
    mutate(position2 = paste(Chr,Start,End,sep="_"))%>% #mutation position without basechange
    mutate(mutID = paste(Sample_orig,position,sep="_")) %>%    #mutation id (here we use Sample_orig, because it is unique)
    mutate(freqID = paste(position,AAChange,sep=":")) %>%   #id to calculate mutation frequency
    mutate(mutFreq = as.vector(table(freqID)[freqID])) %>%  #mutation frequency
    mutate(snp = ((TVAF > 0.4 & TVAF < 0.6)|(TVAF > 0.9))&(AF>0.001))-> tmp
  
  assign(file, tmp)
}

data <- full_join(data1,data2)
#merge variant calling tables from different runs and modify/add columns needed for filtering
df <- data.frame(data) %>% 
  full_join(.,tags) %>% #join with list of tags from manual inspection in igv
  mutate(tag = as.factor(tag)) %>%
  unique %>%
  group_by(Patient.ID,position) %>% 
  mutate(serial.mut = n()) %>% #add column stating whether mutation is present in several samples from the same Patient.ID
  mutate(Genomic_Coordinate_hg38 = paste(Chr,paste("g",Start,sep = "."),paste(Ref,Alt,sep =">"),sep=":"))%>% #Genomic coordinate to match with BRCA_exchange database
  data.frame %>%
  filter(!is.na(Patient.ID)) #filter out Patients without ID 

  as.data.frame(df %>% 
                  group_by(freqID) %>% 
                  summarise(mutID,n(),TVAF,median(TVAF),mean(TVAF),sd(TVAF))) %>% #for every mutation position (as determined by freqID), calculate median VAF, mean VAF and standard deviation. This will be relevant for positions with high mutation frequency
  dplyr::select(mutID,TVAF,'median(TVAF)','sd(TVAF)','mean(TVAF)') %>% 
  dplyr::rename(med.vaf = 'median(TVAF)',sd.vaf = 'sd(TVAF)',mean.vaf = 'mean(TVAF)') %>%
  left_join(df,.) %>% 
  mutate(dev.med = (TVAF-med.vaf)/sd.vaf) %>% #tells us by how many standard deviations the VAF at the current position deviates from the median
  mutate(p.binom = ifelse(mutFreq > 9,log10(1-pbinom(TR2,readDepth,med.vaf)),-Inf))%>%   #p.binom is a measure for the likelihood to get these number of variant reads TR2 (or more), assuming a Bernoulli experiment with p=med.vaf
  filter(!is.na(TR2),!is.na(readDepth))%>% 
  group_by(Patient.ID,position) %>% 
  mutate(serial.mut = n()) %>%
  data.frame %>%
  mutate(hotspot = is.element(position2,hotspots)) %>% #hotspots (defined above)
  mutate(igv = paste(Chr,Start,sep=":"))%>%
  mutate(Patmut=paste(Patient.ID,position,sep="_"))->df
  
rm(data1)
rm(data2)


########   Check for mutations with high frequency but low p.binom--------------------------------
  df %>% filter(!snp,mutFreq>9,p.binom <= -Inf) %>% 
    filter(Func=="exonic",ExonicFunc != "synonymous SNV") %>% 
    filter(med.vaf < 0.44) %>%
    dplyr::select(Sample,cf, Chr, Start, End, Alt, Ref, Gene, ExonicFunc, readDepth, TR1, TR2, TVAF, mutFreq, med.vaf,p.binom,tag,cf) -> df.highfreq 
  
  #filename = "highfreq_variants.xlsx"
  #write.xlsx(df.highfreq,filename,sheetName="somatic",append=TRUE)
  

########   Identify GERMLINE mutations-----------------------------------------------------------
##Determine BRCA status
df %>% 
  filter(Gene == "BRCA1"|Gene =="BRCA2") %>%
  filter(TVAF > 0.25) %>%
  filter(ExonicFunc!= "synonymous SNV") %>%
  filter(AF < 0.05) %>%
  left_join(., brcaexchange, by="Genomic_Coordinate_hg38") %>%
  filter(is.element(BRCA.Exchange_Pathogenicity_expert,c("Pathogenic","Not Yet Reviewed"))|
           (is.na(BRCA.Exchange_Pathogenicity_expert)&
              (is.element(ExonicFunc,c("frameshift substitution","stopgain"))|
                 is.element(Func,c("splicing","exonic;splicing")))))%>% ## all variants that are classified as pathogenic by expert panel or not yet reviewed, or that have no match in BRCA exchange but are truncating
  filter(!str_detect(BRCA.Exchange_Clinical_Significance_ClinVar,"Benign")|is.na(BRCA.Exchange_Clinical_Significance_ClinVar))-> df.brca_germline

ids %>% 
  filter(!str_detect(Sample_orig,"cf")) %>% ###Position zu Sample geändert
  mutate(brca_germline = ifelse(is.element(Patient.ID,df.brca_germline$Patient.ID),1,0))%>%
  dplyr::select(Patient.ID,brca_germline) %>% 
  unique -> id.brca_germline ##this is a list of patient ids with brca status

##identify other HRD Gene germline mutations
hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")

df %>% 
  filter(is.element(Gene,hrd_genes)) %>%
  filter(TVAF > 0.25) %>%
  filter(ExonicFunc!= "synonymous SNV") %>%
  filter(Func == "exonic")%>%
  filter(AF < 0.01) %>%
  filter(is.element(ExonicFunc,c("frameshift substitution","stopgain"))) -> df.hrd_germline#this list has to be discussed with an expert


#Determine IL-6R SNP status
ids %>% filter(!str_detect(Sample_orig,"cf")) %>%
  left_join(.,df %>% 
              filter(Gene == "IL6R") %>% 
              filter(Genomic_Coordinate_hg38 == "chr1:g.154454494:A>C") %>%
              dplyr::select(Sample,Patient.ID,TVAF))%>%
  mutate(il6r_snp = ifelse(TVAF > 0.9,2,
                           ifelse(TVAF < 0.1,0,1))) %>% 
  mutate(il6r_snp = ifelse(is.na(il6r_snp),0,il6r_snp))%>%
  #mutate(Patient.ID = Patient.ID.x)%>%
  dplyr::select(Patient.ID,il6r_snp) %>% 
  unique -> id.il6

########   FILTERING CH calls------------------------------------------------------------
# filter criteria
#functional criteria
df %>%
  filter(ExonicFunc != "synonymous SNV") %>%
  filter(Func == "exonic"|Func == "splicing"|Func == "exonic;splicing") %>%
  dplyr::select(mutID) %>% mutate(mutID = as.character(mutID)) -> mutID.func

## read count criteria
df %>%
  filter(readDepth > 50, TR2 > 9, TVAF > 0.01) %>%
  dplyr::select(mutID) -> mutID.count

## quality criteria
df %>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  dplyr::select(mutID)-> mutID.qual

## mutation call frequency criteria
df %>% 
  filter(AF<0.1) %>%    #filtert alle h?ufig in Datenbanken (=seq errors) gelisteten mutation raus
  filter((mutFreq < 10)|((p.binom<= -10)&med.vaf < 0.44))%>% #filtert alle mutationen raus, die in mehr als 20% der samples auftreten
  dplyr::select(mutID) -> mutID.freq

# rescue ASXL1 dupG mutations <- this step is no longer needed, when we use p.binom 
#df %>%
#  filter(str_detect(Start,'32434638'))%>%filter(str_detect(AAChange,'ASXL1')) %>%
#  mutate(dev.med = ((TVAF - median(TVAF))/sd(TVAF))) %>% #calculate deviation from median in terms of standarddeviations
#  filter(dev.med > 1) %>%    #to be discussed
#  dplyr::select(mutID) -> mutID.asxl1

## rescue hotspots
df %>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 3) %>%
  filter(TVAF >0.005) %>%
  filter(is.element(position,hotspots))%>%
  dplyr::select(mutID)-> mutID.hotspots

## rescue known CHIP mutations with weakened quality criteria
df %>% filter(ChipPub != "") -> mutID.CHIP
df %>%
  filter(FisherScore < 20) %>% 
  filter(StrandBalance2 != 1 & StrandBalance2 != 0) %>%     #filter out mutations only seen on one strand
  filter(TR2 > 7) %>%
  filter(TVAF >0.005) %>%
  dplyr::select(mutID)-> mutID.CHIP.qual


#filtering
##somatic variants
inner_join(mutID.func,mutID.count) %>% 
  inner_join(.,mutID.freq) %>% 
  inner_join(.,mutID.qual) %>% 
  inner_join(.,df) %>% 
  full_join(.,inner_join(df,inner_join(mutID.CHIP,mutID.CHIP.qual))) %>%
  #full_join(.,inner_join(df,mutID.asxl1)) %>%
  full_join(.,inner_join(df,mutID.hotspots))%>%
  filter(snp == FALSE) %>%
  mutate(current_filter = 1) %>% ##tag all variants passing current filter, then join with list of previously tagged true
  full_join(df %>% filter(tag=="true"))-> df.filtered

#rescue Step (when Mutation is found in more than one timepoints
df.filtered%>%
  filter(tag== "true")%>%
  select(.,Patmut)%>%
  unique()->Mutations
  Mutations<-unlist(Mutations)
for (i in Mutations) {
  print(i)
  df%>%
    filter(Patmut == i)->x
  nrow(x)->y
  if (y >= 2)
  {df%>%
      filter(Patmut == i)->x
      x%>%
        mutate(tag = "true")->x
    bind_rows(df.filtered, x)%>%
      unique()->df.filtered}
  if (y == 1)
  {}
  rm(x)
  rm(y)
}
rm(i)
rm(Mutations)

df.filtered %>% filter(Visite == "C1D1") -> df.filtered.c1d1


#write into excel file

#filename = "filtered_results.xlsx"
#write.xlsx(df.filtered,filename,sheetName="test",append=TRUE)
#write.xlsx(df.brca_germline,filename,sheetName="BRCA",append=TRUE)
#write.xlsx(df.hrd_germline,filename,sheetName="HRD",append=TRUE)
#write.xlsx(df.syn,filename,sheetName="Synonymous",append=TRUE)
#write.xlsx(df,filename,sheetName="all",append=TRUE)

##SNP
inner_join(mutID.count,mutID.qual) %>% 
  inner_join(.,df) %>% 
  filter(snp == TRUE)%>%
  filter(avsnp150!=".")%>%
  filter(FisherScore>1)%>%
  filter(AF >= 0.3)-> df.snp

rm(mutID.CHIP)+
rm(mutID.CHIP.qual)+
rm(mutID.count)+
rm(mutID.freq)+
rm(mutID.func)+
rm(mutID.hotspots)+
rm(mutID.qual)+
rm(tags)+
rm(tmp)+
rm(file)+
rm(data)+
rm(brcaexchange)
########   SNP-Comparison -> Heatmap ####

QC_germline<-read.table(file = 'QC/Total_result.txt', sep = '\t', header = FALSE)
names(QC_germline) <- c("Sample1","Sample2","Score", "MatchingStatus")
filename="QC/QC_germline.xlsx"
#write.xlsx(QC_germline,filename, sheetName="QC", append=TRUE)

##Preparation for Heatmap
df.snp%>%
  mutate(PatSample = paste(Patient.ID,Sample,Visite,sep="_"))->df.snp
##rescue matching germline Mutations
select(df.snp,Patmut)->b
b=unlist(b)
for (Mut in b)
{print(Mut)
  df%>%
    filter(Patmut == Mut)->a
  bind_rows(df.snp,a)->df.snp
  unique(df.snp)->df.snp}
#oder?
select(df.snp,Patmut)->q
semi_join(df, df.snp, by="Patmut")->df.snp

q->df.snp
rm(q)
#speichern als backup
df.snp->df.snpbackup

df.snp%>%
  filter(Patient.ID != "0")->df.snp
#Heatmap (Robert Skript)
### long table to wide table
newdata <- dcast(data = df.snp,
                 formula = position ~ PatSample,    
                 length)         

### make first column the rowname
rownames(newdata) <- newdata[, "position"]
### and remove it
newdata <- newdata[,-1]
#colSums(newdata)->a

### create distance matrix (type "manhatten")
dst <- dist(t(newdata), method = "manhattan")
dst <- data.matrix(dst)

### plotting
dim <- ncol(dst)
pdf(file="QC/Distanz.pdf",height=35,width=45)
image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, rownames(dst), cex.axis = 0.5, las=3)
axis(2, 1:dim, rownames(dst), cex.axis = 0.5, las=1)

text(expand.grid(1:dim, 1:dim), sprintf("%1i", dst) , cex=0.6)
dev.off()

#write newdata in excelfile
filename="QC/Distanz.xlsx"
write.xlsx(newdata,filename,sheetName=" ",append=TRUE)

rm(newdata)
rm(dim)
rm(dst)

##comparison germline:
# Get the unique patient IDs
patient_ids <- unique(df.snp$Patient.ID)

# Initialize a data frame to store the results
germline_comparison <- data.frame(patient_id = character(0), sample_id_1 = character(0), sample_id_2 = character(0), overlap_percentage = numeric(0))

# Loop through each patient ID
for (patient_id in patient_ids) {
  
  # Subset the data frame to only include samples for the current patient
  patient_samples <- df.snp[df.snp$Patient.ID == patient_id, ]
  
  # Get the unique sample IDs for the patient
  sample_ids <- unique(patient_samples$Sample.ID)
  
  # If there is only one sample for the patient, skip the patient
  if (length(sample_ids) < 2) {
    next
  }
  
  # Loop through each sample combination
  for (i in 1:(length(sample_ids)-1)) {
    for (j in (i+1):length(sample_ids)) {
      
      # Subset the data frame to only include the current samples
      sample_1 <- patient_samples[patient_samples$Sample.ID== sample_ids[i], ]
      sample_2 <- patient_samples[patient_samples$Sample.ID == sample_ids[j], ]
      
      # Get the intersection of the two samples
      intersection <- intersect(sample_1$position, sample_2$position)
      
      # Calculate the overlap percentage
      overlap_percentage <- (length(intersection) / nrow(sample_1)) * 100
      
      # Store the results in the results data frame
      germline_comparison <- rbind(germline_comparison, data.frame(patient_id = patient_id, sample_id_1 = sample_ids[i], sample_id_2 = sample_ids[j], overlap_percentage = overlap_percentage))
    }
  }
}

# View the results
germline_comparison

rm(patient_id)+
rm(patient_ids)+
rm(patient_samples)+
rm(sample_1)+
rm(sample_2)+
rm(i)+
rm(j)+
rm(intersection)+
rm(overlap_percentage)+
rm(sample_ids)
  

########   Mutational Analysis------------------------------------------------------------------------------

##PLOTs

#Gene Mutation Prevalence Plot (plots number of gene-x-mutated patients)
ids%>%
  filter(Patient.ID!=0)%>%
  select(.,Patient.ID)%>%
  unique()->ID
nrow(ID)->nop
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
  theme_gray() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  coord_flip() + 
  scale_fill_manual(values = c("non HRD" = "#486081", "HRD" = "#88acd4")) -> p.mutprev
p.mutprev

png("plots/mutprev.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev
dev.off()

#CH prevalence by BRCA status
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  dplyr::select(Patient.ID) %>% unique -> id.ch
id.brca_germline %>% 
  mutate(CH = ifelse(is.element(Patient.ID,id.ch$Patient.ID),1,0))%>%
  dplyr::select(brca_germline,CH) %>% 
  table 


##prevalence plot by BRCA status
#prevalences in BRCA wildtype patients
nop <- 87-19
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
  mutate(prev = Freq/nop) %>% 
  arrange(prev) %>%
  mutate(brca = 0) -> prev.table_brca0
names(prev.table_brca0)<- c( "Gene","Freq","prev","brca")

#prevalences in brca mutated patients
nop <- 19 #number of brca mut patients
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
  mutate(prev = Freq/nop) %>% 
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
  theme_graphicalabstract() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  scale_fill_viridis(discrete=TRUE) +
  scale_fill_manual(values = c("0" = "#486081", "1" = "#88acd4")) +
  coord_flip() -> p.mutprev
p.mutprev
png("plots/mutprev-BRCA.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev
dev.off()
########   Create txt.file for Circleplot ####

###df nur mit true und Tag 1 Ergebnissen
df.filtered%>%
  filter(Visite == "C1D1")%>%
  filter(tag== "true")%>%
  filter(Patient.ID !=0)%>%
  filter(cf == "0")->data.frame

##Gene rausfiltern
select(data.frame,Gene)->genes
unique(genes)->genes ->genes2
genes<-unlist(genes)

#create tables to fill
genes2->finaltable

tab <- matrix(c(0,0),ncol=2, byrow=TRUE)
colnames(tab) <- c('y','X')
tab <- as.table(tab)

#Create table
for (b in genes)
{(print(b))
  data.frame%>%
    filter(Gene == b)->df.gene
  nrow(df.gene)->f
  select(df.gene,Patient.ID)->geneID
  inner_join(data.frame,geneID, by="Patient.ID")->df.geneID
  df.geneID%>%
    filter(Gene != b)->df.geneID
  genes2%>%
    filter(Gene != b)->genesx
  genesx=unlist(genesx)
  
  for (c in genesx)
  {df.geneID%>%
      filter(Gene == c)->a
    nrow(a)->a
    
    tab2 <- matrix(c(c,a),ncol=2, byrow=TRUE)
    tab2 <- as.table(tab2)
    
    colnames(tab) <- c(b,'X')
    rbind(tab,tab2)->tab}
  
  as_tibble(tab)->Col1
  colnames(Col1) <- c("Gene", b)
  
  left_join(finaltable,Col1,by="Gene")->finaltable
  
  tab <- matrix(c(0,0),ncol=2, byrow=TRUE)
  colnames(tab) <- c('y','X')
  tab <- as.table(tab)
  
  #doppelte Zeilen entfernen  
  finaltable<-unique(finaltable)
}

##Sums an seite schreiben
#create Matrix
tab3<-matrix(c(0,0),ncol=2,byrow=TRUE)
colnames(tab3)<-c('data','Gene')
tab3<-as.table(tab3)

for (b in genes)
{data.frame%>%
    filter(Gene == b)->df.gene
  nrow(df.gene)->f
  
  tab4<-matrix(c(f,b),ncol=2,byrow=TRUE)
  colnames(tab3)<-c('data','Gene')
  tab4<-as.table(tab4)
  
  rbind(tab3,tab4)->tab3
}

as_tibble(tab3)->tab3

tab3%>%
  filter(data != "0")->tab3

##join Gene Sums mit finaltable
left_join(tab3,finaltable, by="Gene")->finaltable


### for creating first two columns
finaltable->backup
select(backup,data,Gene)->newROW

##create and adapt first two columns
matrix(c("data","data","data","Gene"),ncol=2,byrow=TRUE)->newROW2
as.matrix(newROW)->newROW
rbind(newROW2,newROW)->newCOL

##tilt
as.data.frame(t(finaltable))->finaltable

##Bind to final table
bind_cols(newCOL,finaltable)->finaltable

###replace NA by 0
finaltable[is.na(finaltable)] <- 0


#Wirte into excel file
filename="Circleplot.xlsx"
write.xlsx(finaltable,filename,sheetName="X",append=TRUE)

#Write into txt file
write.table(finaltable, file = "plots/Circleplot.txt", sep = " ",
            row.names = FALSE, col.names = FALSE)

rm(backup)+
rm(genesx)+
rm(tab)+
rm(tab3)+
rm(tab4)+
rm(finaltable)+
rm(genes)+
rm(f)+
rm(c)+
rm(b)+
rm(a)+
rm(newROW)+
rm(newROW2)+
rm(newCOL)+
rm(df.geneID)+
rm(df.gene)+
rm(Col1)+
genes2->genes+
rm(genes2)

########   GenVisar - waterfall plot ####
library(viridisLite)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenVisR")

library(GenVisR)


select(data.frame,Patient.ID,Gene,ExonicFunc)->Target_Seq
colnames(Target_Seq) <- c("Patient.ID","Gene","Consequence")

## (.) mit splicesite erstetzen
select(Target_Seq,Consequence)->test2
select(Target_Seq,Patient.ID,Gene)->test3

test2<-replace(test2$Consequence,test2$Consequence==".", "splicesite")
test2<-data.frame(test2)
colnames(test2) <- "Consequence"

bind_cols(test3,test2)->Target_Seq

# List of Genes in panel (will be needed for the gene.order vector)
# load genes
colnames(genes) <- "Hugo_Symbol"

# Select columns for landscape

# create MAF file as input for waterfall() function
# needs to be in data frame format 

OvCAmaf <- Target_Seq$Patient.ID
OvCAmaf <- as.data.frame(OvCAmaf)                                         
names(OvCAmaf) <- "Patient.ID"
OvCAmaf$Patient.ID <- as.character(OvCAmaf$Patient.ID)

OvCAmaf$Hugo_Symbol <- Target_Seq$Gene
OvCAmaf$Variant_Classification <- Target_Seq$Consequence 

# count sample number
n <- length(unique(OvCAmaf$Patient.ID)
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
#data.frame(OvCAmaf)->OvCAmaf
waterfall(OvCAmaf, 
          fileType = "Custom",            ##fileType can be "MAF", "MGI", "Custom"
          rmvSilent = FALSE,               ## don't show silent mutations (we don't have any)
          mainDropMut = TRUE,             ## unused mutation types will be dropped from the legend
          mainPalette = viridis(5),          ## previously specified colors
          plotMutBurden = FALSE,          ## mutation burden plot at the top
          mainLayer = mainLayer,          ## specify defined layers
          sampRecurLayer = sampRecurLayer, 
          geneOrder = gene.order,         ## order to plot the genes (default: decreasing gene freq.)
          variant_class_order = variants, 
          section_heights = c(0.15, 1),   ## relative size of the different plots
          mainXlabel = FALSE,
          mainLabelSize = 0,
          main_geneLabSize = 1)->p.waterfall
##ich kann das nicht als png speichern?!
png("plots/waterfall.png",width=6, height=6,units="in",res=500,type="cairo")
p.waterfall
dev.off

rm(variants)+
rm(tab2)+
rm(p.waterfall)+
rm(gene.order)+
rm(test3)+
rm(test2)+
rm(Target_Seq)+
rm(sampRecurLayer)+
rm(OvCAmaf)+
rm(mainLayer)+
rm(genes)

########   number of mutations plot------------------------------------------------------------
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%  
  dplyr::select(Gene,ExonicFunc) %>% 
  mutate(ExonicFunc = replace(ExonicFunc,ExonicFunc == ".","splice mutation")) %>%
  data.frame %>% 
  group_by(Gene) %>% 
  summarise(Gene.freq = n())%>%
  as.data.frame %>% 
  full_join(df.filtered %>% filter(tag == "true"),.,by = "Gene") -> df.genefreq

##mutation frequency according to type of mutation (not yet working properly)
df.genefreq %>%
  dplyr::select(Gene,ExonicFunc) %>% 
  mutate(ExonicFunc = replace(ExonicFunc,ExonicFunc == ".","splice mutation")) %>%
  data.frame %>% 
  #select(Gene) %>% 
  table %>% 
  data.frame %>% 
  filter(Freq >0) %>% 
  inner_join(.,df.genefreq %>% 
               dplyr::select(Gene,Gene.freq) %>% 
               unique(),by = "Gene") %>%
  arrange("Gene.freq")  %>% 
  ggplot(aes(x=reorder(Gene, Gene.freq), y=Freq,fill=ExonicFunc)) +
  geom_bar(stat="identity", width=0.6)+
  #geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(position = "right")+
  ylab("Mutation Frequency") +
  theme_graphicalabstract() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  scale_fill_viridis(discrete=TRUE) +
  labs(fill="")+
  coord_flip() -> p.mutfreq
p.mutfreq

ggsave("H:/Meine Ablage/plots/mutfreq.png")
#oder
png("plots/mutfreq.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutfreq
dev.off()

rm(p.mutfreq)+
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
