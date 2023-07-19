# ==============================================================================
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: Creates a MAF-like format from our inhouse pipeline output that can be used with MAFtools etc
#
# Input: finalized variant calling table as data frame containing:
#         Func, AAChange,
#
# Output: data.frame in MAF-like format
#
# ==============================================================================
########   Dependencies   #####
library(base)
library(dplyr)
library(xlsx)
library(stringr)
library(reshape)
library(tidyr)
library(readxl)
library(reshape2)
library(forcats)

########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

makeMAF <- function(df){
df %>% 
  filter(Func == "exonic"|Func=="splicing")%>%
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
  dplyr::select(Gene,amino_acid_change,AA_position,Patient.ID,ExonicFunc,Chr,Start,End,Ref,Alt,TVAF)->df.temp

names(df.temp) <- c(
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
  "VAF"
)

variant_replacements = c('nonframeshift substitutiondel' = 'In_Frame_Del',
                      'nonframeshift substitutionins' = 'In_Frame_Ins',
                      'nonframeshift substitution.' = 'Missense_Mutation',
                      'nonsynonymous SNV.' =  'Missense_Mutation',
                      'stopgain.' = 'Nonsense_Mutation',
                      'stopgaindel' =  'Nonsense_Mutation',
                      'stopgainins' = 'Nonsense_Mutation',
                      'stoploss.' = 'Nonstop_Mutation',
                      'stoplossdel' = 'Nonstop_Mutation',
                      'frameshift substitutiondel' = 'Frame_Shift_Del',
                      'frameshift substitutionins' =  'Frame_Shift_Ins',
                      '..' = "Splice_Site")

indel_replacements = c('ins' = 'INS', 
                       'del' = 'DEL',
                       '.' = 'SNP')

df.temp %>% mutate(Tumor_Seq_Allele2 = Variant_Allele,
                    Variant_Classification = Mutation_Type,
                    Variant_Type = Mutation_Type,
                    Tumor_Sample_Barcode = Sample_ID,
                    NCBI_Build = "GRCh38") %>%
  mutate(indel = ifelse(nchar(Variant_Allele)>nchar(Reference_Allele),"ins",
                        ifelse(nchar(Variant_Allele)<nchar(Reference_Allele),"del",".")))%>%
  mutate(Variant_Type=paste(Mutation_Type,indel,sep="")) %>%
  mutate('Variant_Classification' = variant_replacements[Variant_Type]) %>%
  mutate('Variant_Type' = indel_replacements[indel])%>%
  dplyr::select(-indel)-> output
print(output)
}