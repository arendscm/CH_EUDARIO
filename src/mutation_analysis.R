# ______________________________________________________________________________
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: Mutation analysis and plots from filtered mutational data
#
# Input: df.filtered from data/interim/seqdata_filtered.RData
#
# Output: Excel list of filtered results, plots, ...
# press ALT-O
# ______________________________________________________________________________
########   Dependencies and load data  ####
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
library(RColorBrewer)
library(ggsci)

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata_filtered.RData')
df.filtered.c1d1 <-df.filtered.c1d1[is.na(df.filtered.c1d1$replicate),]

######## Get Patient ids
source("src/material_table.R")

######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")
source("src/genegroup_definitions.R")



######## Available samples Overview  ##############

df.material %>% 
  dplyr::select(-sum_wb,-sum_cf)%>%
  mutate(eot_ue_cf=ue_cf+eot_cf)%>%
  dplyr::select(-eot_cf,-ue_cf)%>%
  melt.data.frame(id="Patient.ID") %>%
  ggplot(aes(y=variable,x=Patient.ID,fill=as.factor(value)))+
  geom_tile()+
  scale_fill_npg(name="Sample sequenced")+
  xlab("Visit")+
  ylab("Patient ID")+
  my_theme()+
  theme(axis.text.x=element_text(angle=90))->p.material

png("output/figures/material.png",width=15, height=7,units="in",res=500,type="cairo")
p.material
dev.off()

#get numbers
df.material %>% mutate(eot_ue_cf = eot_cf+ue_cf) %>% dplyr::select(-Patient.ID) %>% sign() %>% colSums()

########   Mutational Analysis preparation #####

##PLOTs
df.filtered.c1d1 %>% mutate(Gene=ifelse(Gene=="U2AF1;U2AF1L5","U2AF1",Gene)) -> df.filtered.c1d1

nop <- ids%>%
  filter(Visite == "C1D1" & Material == "wb")%>%
  select(.,Patient.ID)%>%
  unique()%>%nrow #number of patients

#number of CH positive Patients
df.filtered.c1d1%>%
  filter(tag == "true" & TVAF >= 0.01)%>%
  select(.,Patient.ID)%>%
  unique()%>%
  nrow()

#number of mutations
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>% 
  nrow() #number of mutations


########   Gene Mutation Prevalence Plot (plots number of gene-x-mutated patients)  #####
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
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
  mutate(DDR = ifelse(is.element(Gene,ddr_genes),"DDR","non-DDR"))%>%
  ggplot(aes(x=reorder(Gene, Freq), y=prev, fill=DDR)) +
  geom_bar(stat="identity", width=0.6)+
  geom_text(aes(label=Freq), hjust= -1, vjust=0.35, size=4)+
  xlab("")+
  scale_y_continuous(labels = percent,limits=c(0,0.35), position = "right")+
  ylab("Gene Mutation Prevalence [%]") +
  my_theme() +
  theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.35,face="italic")) +
  coord_flip() + 
  scale_fill_npg() -> p.mutprev
p.mutprev

png("output/figures/mutprev.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev
dev.off()


########   Number of Mutations per Gene with type of mutation  ####



#########    No of mutations per patient
ids %>% 
  filter(visit_material=="C1D1_wb") %>% 
  dplyr::select(Patient.ID) %>% 
  left_join(.,df.filtered.c1d1%>%
              filter(TVAF > 0.01,tag=="true")%>%
              group_by(Patient.ID)%>%
              mutate(n.patient = n())%>%
              data.frame%>%
              dplyr::select(Patient.ID,n.patient)%>% 
              unique)%>%
  mutate(n.patient=replace_na(n.patient,0))%>%
  mutate(n.patient=ifelse(n.patient > 5, ">5",n.patient))%>%
  mutate(n.patient=factor(n.patient,levels=c("0","1","2","3","4",">5")))%>%
  dplyr::select(n.patient)%>%
  table%>%
  data.frame -> df.nom
names(df.nom) <- c("nom","Freq")

df.nom %>% 
  filter(nom!=0)%>%
  ggplot(., aes(x = nom, y = Freq, fill="1")) +
  geom_bar(stat = "identity") +
  labs(x = "Number of Mutations", y = "Number of Patients")+
  scale_fill_npg(name="",breaks=c(""))+
  #scale_x_continuous(limits = c(0.5, 8.5), breaks = 1:8)+
  #ggtitle("No. of mutations per patient [>1%]")+
  my_theme()->p.nom
p.nom

png("output/figures/nom.png",width=4, height=3,units="in",res=500,type="cairo")
p.nom
dev.off()


########   mutation frequency according to type of mutation  ####
df.filtered.c1d1 %>% 
  filter(tag == "true") %>%  
  filter(TVAF >= 0.01) %>%
  dplyr::select(Gene,ExonicFunc) %>% 
  mutate(ExonicFunc = replace(ExonicFunc,ExonicFunc == ".","splice mutation")) %>%
  data.frame %>% 
  group_by(Gene) %>% 
  summarise(Gene.freq = n())%>%
  as.data.frame %>% 
  full_join(df.filtered.c1d1 %>% filter(tag == "true"),.,by = "Gene") -> df.genefreq

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
  scale_fill_npg() +
  labs(fill="")+
  coord_flip()+
  theme(legend.position = c(0.7, 0.3))-> p.mutfreq
p.mutfreq


png("output/figures/mutfreq.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutfreq
dev.off()



########   MAFTOOLS - Oncoplot #### 

library(maftools)

###turn df.filtered into MAF compatible format
df.maf <- read.maf(makeMAF(df.filtered.c1d1%>% filter(tag=="true",TVAF >= 0.01)))

oncoplot(df.maf)


########   Rose Chart ####
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
  
  annotate("text", x = rep(max(data$id),4), y = c(0,0.05,0.1,0.2), label = c("0","5%","10%","20%") , color="black", size=2 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(stat="identity", width=0.9) +
  
  ylim(-0.1,0.35) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank())+
  coord_polar()+
  
  geom_segment(data=base_data, aes(x = start, y = -0.01, xend = end, yend = -0.01), colour = "grey", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = -0.032, label=Gene), colour = "black", alpha=1, size=2, fontface="bold", inherit.aes = FALSE)


p.rosechart

png("output/figures/p.rosechart.png",width=8, height=8,units="in",res=500,type="cairo")
p.rosechart
dev.off()



########   prevalence plot by PATHOGENIC BRCA status####
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
load('data/interim/id.BRCA_path.RDATA')
#How many BRCA+/- and CH+/-
id.brca_germline_path%>%
  filter(brca_germline == 1)%>%
  .$Patient.ID->ID.BRCA.path  # n() -> how many are BRCA mutated
#not BRCA Muatated patients = 94-ID.BRCA.path
df.filtered.c1d1%>%
  filter(Patient.ID %in% ID.BRCA.path)%>%
  filter(TVAF >= 0.01)%>%
  filter(tag == "true")->mutations_in_BRCA_pos
mutations_in_BRCA_pos$Patient.ID%>%
  unique->No.ID.pos #how many are BRCA+ and CH+
df.filtered.c1d1%>%
  filter(!Patient.ID %in% ID.BRCA.path)%>%
  filter(TVAF >= 0.01)%>%
  filter(tag == "true")->mutations_in_BRCA_neg
mutations_in_BRCA_neg$Patient.ID%>%
  unique->No.ID.pos #how many are BRCA- and CH+
n.brcamut <- id.brca_germline_path %>% 
  mutate(CH = ifelse(is.element(Patient.ID,id.ch$Patient.ID),1,0))%>%
  dplyr::select(brca_germline) %>% sum 

df.filtered.c1d1 %>% 
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  left_join(.,id.brca_germline_path,by = "Patient.ID")%>%
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
  left_join(.,id.brca_germline_path,by = "Patient.ID")%>%
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
  coord_flip() -> p.mutprev_path
p.mutprev_path

png("output/figures/mutprev-BRCA_path.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev_path
dev.off()



########   multiple mutations plots  ####
df.filtered.c1d1%>%
  filter(tag == "true" & TVAF >= 0.01)->test

test<-test%>%
  group_by(Patient.ID) %>%
  filter(n() >= 3) %>%
  ungroup()

Comutations_all <-ggplot(test, aes(x = position, y = TVAF, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Patient.ID, scales = "free_x") +
  labs(x= "",y = "TVAF", fill = "Genes")+
  theme(axis.text.x = element_blank(),
        xlab = NULL,
        axis.text.y = element_text(size = 8))+
  scale_y_continuous(breaks = c(0.01, 0.1, 0.2, 0.3),
                     labels = c(0.01, 0.1, 0.2, 0.3))
#scale_fill_brewer(palette = "Set3")
png("output/figures/comutationover3.png",width=6, height=6,units="in",res=500,type="cairo")
Comutations_all
dev.off()


########   PPM1D Comutations:  #####
df.filtered.c1d1%>%
  filter(tag == "true" & TVAF >= 0.01)->test

filtered_df <- test %>%
  filter(Gene == "PPM1D")

# Filter the data to include only patients with at least two mutations in the gene "PPM1D"
filtered_df <- filtered_df %>%
  group_by(Patient.ID) %>%
  filter(n() >= 2) %>%
  ungroup()

filtered_df$Patient.ID%>%
  unique()->Patient.ID.PPM1DCo

df.filtered.c1d1%>%
  filter(tag == "true" & TVAF >= 0.01)->test

test<-test%>%
  group_by(Patient.ID) %>%
  filter(n() >= 2) %>%
  ungroup()

test%>%
  filter(Patient.ID %in% Patient.ID.PPM1DCo)->test

p.PPM1DCo<- ggplot(test, aes(x = position, y = TVAF, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Patient.ID, scales = "free_x") +
  labs(x= "",y = "TVAF", fill = "Genes")+
  theme(axis.text.x = element_blank(),
        xlab = NULL,
        axis.text.y = element_text(size = 8))+
  scale_y_continuous(breaks = c(0.01, 0.1, 0.2, 0.3),
                     labels = c(0.01, 0.1, 0.2, 0.3))
#scale_fill_brewer(palette = "Set3")

png("output/figures/PPM1D-Comutations.png",width=6, height=6,units="in",res=500,type="cairo")
p.PPM1DCo
dev.off()
