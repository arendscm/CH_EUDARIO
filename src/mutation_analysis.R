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

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata_filtered.RData')
df.filtered.c1d1 <-df.filtered.c1d1[is.na(df.filtered.c1d1$replicate),]

######## Get Patient ids
source("src/ids.R")
ids <- ids[is.na(ids$replicate), ]


######## Functions and themes
source("src/createMAF.R")
source("src/global_functions_themes.R")

########   Mutational Analysis preparation #####
##PLOTs
nop <- ids%>%
  filter(Visite == "C1D1" & Material == "wb")%>%
  select(.,Patient.ID)%>%
  unique()%>%nrow #number of patients

df.filtered.c1d1%>%
  filter(tag == "true" & TVAF >= 0.01)%>%
  select(.,Patient.ID)%>%
  unique()%>%
  nrow()->no.chpos.pat


hrd_genes <- c("ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3", "BRCA1", "BRCA2")
ch_genes <- c("DNMT3A", "TET2" ,  "JAK2" ,  "ASXL1" , "SF3B1" , "SRSF2" , "TP53"  , "U2AF1" , "PPM1D" , "CBL"  ,  "IDH1"  , "IDH2"  , "BCOR"  , "BCORL1", "EZH2" ,  "RAD21" , "STAG2" , "CHEK2" , "GNAS"  , "GNB1"  , "ATM"   , "KRAS" ,  "NRAS",   "WT1" ,   "MYD88" ,
              "STAT3" , "BRCC3" , "CALR"  , "CEBPA" , "CSF3R" , "ETV6"  , "FLT3" ,  "GATA2" , "GATA1" , "KIT" ,   "MPL" ,   "NPM1" ,  "PTPN11" ,"RUNX1" , "SETBP1" ,"NF1"  ,  "PHF6")

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

png("output/figures/mutprev.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutprev
dev.off()

rm(prev.table)
rm(p.mutprev)

########   Number of Mutations per Gene   ####
# create a summary data frame showing the number of mutations per gene
mutation_counts <- df.filtered.c1d1%>%
  filter(tag == "true" & TVAF >= 0.01)%>%
  group_by(Gene) %>%
  summarize(count = n())%>%
  arrange(desc(count))

#create a barplot of the mutation counts per gene
ggplot(mutation_counts, aes(x = Gene, y = count, fill = "Mutation Counts")) +
  geom_bar(stat = "identity", fill = "#486081") +
  geom_text(aes(label = count), vjust = -0.5) + # add labels for mutation counts
  #coord_flip() + # flip the coordinates to create a horizontal barplot
  xlab("") +
  ylab("Number of mutations") +
  ggtitle("Mutation counts per gene") +
  scale_x_discrete(limits = mutation_counts$Gene) + # set x-axis limits based on sorted order
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p.no.mut

png("output/figures/no.mut.png",width=6, height=6,units="in",res=500,type="cairo")
p.no.mut
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

########   No of mutations per patient ####
df.filtered.c1d1%>%
  filter(tag == "true" & TVAF >= 0.01)->test
mut_count <- test %>%
  group_by(Patient.ID) %>%
  summarize(mutations = n())

mutation_barplot <- mut_count %>%
  group_by(mutations) %>%
  summarize(num_patients = n())

# Create the bar plot with ggplot2
ggplot(mutation_barplot, aes(x = mutations, y = num_patients)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Number of Mutations", y = "Number of Patients")+
  scale_x_continuous(limits = c(0.5, 8.5), breaks = 1:8)+
  ggtitle("No. of mutations per patient [>1%]")->p.nom
p.nom

png("output/figures/no of mutations.png",width=6, height=6,units="in",res=500,type="cairo")
p.nom
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

########   GenVisar - waterfall plot #### 

library(viridisLite)

library(GenVisR)

###turn df.filtered into MAF compatible format
df.maf <- makeMAF(df.filtered.c1d1%>% filter(tag=="true"))
select(df.maf,Hugo_Symbol)%>%unique()->gene
#save(gene,file="data/interim/gene_Waterfall.RDATA")
#or
  load('data/interim/gene_Waterfall.RDATA')

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

gene.order <- gene$Hugo_Symbol                      ## vector specifying target genes and gene order

variants <- unique(df.maf$Variant_Classification)    ## vector specifying types of mutation

# or
# variants <- c("Nonsynonymous SNV",                 ## you can choose any type of mutation you want
#              "Frameshift Deletion",                ## for the "costum" file type
#             "Splice Site", 
#              "Stop-Gain"
#               )

# Generate plot and save as .pdf
#cant save the plot...it will just show white

########   number of mutations plot  (nom) ------------------------------------------------------------
df.filtered.c1d1 %>%
  filter(tag == "true") %>%
  filter(TVAF >= 0.01) %>%
  dplyr::select(nom)%>% 
  table %>% 
  data.frame %>%
  filter(. != 0)%>%
  ggplot(aes(y=Freq,x=.,fill=.)) +
  geom_bar(stat="identity", width=0.6)+
  ylab("Number of Patients") +
  xlab("Number of Mutations")+
  scale_fill_manual(values=c("#31688EFF","#21908CFF","#21908CFF","#21908CFF"),name="",labels=NULL)+
  theme_Publication() +theme(legend.position = "none")-> p.nom
p.nom

png("nom05.png",width=4, height=4,units="in",res=500,type="cairo")
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
  scale_fill_viridis(discrete=TRUE) +
  labs(fill="")+
  coord_flip() -> p.mutfreq
p.mutfreq


png("output/figures/mutfreq.png",width=6, height=6,units="in",res=500,type="cairo")
p.mutfreq
dev.off()

rm(p.mutfreq)
rm(df.genefreq)

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

