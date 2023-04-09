# ______________________________________________________________________________
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: determine mutation signatures in cfDNA mutations and wb mutations. 
#
# Input: df from data/interim/seqdata.RData, ...
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
library(rtracklayer)
library(GenomicRanges)
library(liftOver) #tool to liftover between genomic assemblies (hg19, hg38)
library(deconstructSigs) #tool to extract trinucleotide context, and decompose contributation to mutation signature
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
#####  Data preparation ####
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
Categories<-c('CH','HRD','other','TP53')

#dataframe with mutation calls from cfDNA             
df %>% 
  filter(Material=="cf") %>% 
  filter(Visite == "C1D1")%>%
  filter(is.na(replicate))%>%
  filter(!is.element(Sample.ID, failedSamples))%>%
  dplyr::select(all_of(variables)) %>%
  mutate(cfID=paste(Patient.ID,position,sep="_"))-> df.cf

#data frame with mutation calls from WB samples that have matched cfDNA samples
df %>% 
  filter(is.element(Patient.ID,df.cf$Patient.ID)) %>% 
  filter(is.na(replicate))%>%
  filter(Material=="wb") %>% 
  filter(Visite == "C1D1")%>%
  filter(!is.element(Sample.ID, failedSamples))%>%
  dplyr::select(all_of(variables)) %>%
  mutate(cfID=paste(Patient.ID,position,sep="_")) -> df.cf_wb

#####  Determine trinucleotide context of mutations cf vs wb



##create data frame SNVs with compartment (cf or wb)

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
  mutate(compartment = ifelse(TVAF.x > TVAF.y*5,"cf","wb"))%>%
  mutate(chr = Chr.x,
         sample.id = Patient.ID.x,
         start = Start.x,
         end = End.x,
         REF = Ref.x,
         ALT= Alt.x,
         type = ExonicFunc.x)%>%
  dplyr::select(c("chr",
                  "start",
                  "end",
                  "REF",
                  "ALT",
                  "compartment"))-> df.sig

df.sig %>%  
  filter(type=="nonsynonymous SNV"|type=="synonymous SNV")->df.sbs


df.sig %>%  
  filter(str_length(REF)==2&str_length(ALT)==2)->df.dbs

  
##create GRanges object for sbs
gr.sig <- makeGRangesFromDataFrame(df.sig,
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)

genome(gr.sig) <- "hg38"


###granges object for dbs 
gr.dbs <- get_mut_type(gr.sig,type="dbs",predefined_dbs_mbs = TRUE)
gr.snv <- get_mut_type(gr.sig,type="snv",predefined_dbs_mbs = TRUE)
gr.indel <- get_mut_type(gr.sig,type="indel",predefined_dbs_mbs = TRUE)

gr.dbs$REF <- DNAStringSet(gr.dbs$REF)
gr.dbs$ALT <- DNAStringSetList(DNAStringSetList(as.list(gr.dbs$ALT)))


##SNV contexts
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
hg19 <- "BSgenome.Hsapiens.UCSC.hg19"


split(gr.snv,as.factor(gr.snv$compartment)) -> grl.comp.sbs
mut_mat.sbs <- mut_matrix(grl.comp.sbs, ref_genome = ref_genome)
plot_96_profile(mut_mat.sbs)-> p

png("output/figures/mut_sbs.png",width=10, height=4,units="in",res=500,type="cairo")
p
dev.off()

##dbs contexts (funktioniert nicht)

split(gr.dbs,as.factor(gr.dbs$compartment)) -> grl.comp.dbs
dbs_context <- get_dbs_context(grl.comp.dbs)
dbs_counts <- count_dbs_contexts(dbs_context)
plot_dbs_contexts(dbs_counts, same_y = TRUE) -> p.dbs

png("output/figures/mut_dbs.png",width=10, height=4,units="in",res=500,type="cairo")
p.dbs
dev.off()



##mut type occurences
type_occurrences <- mut_type_occurrences(grl.comp.sbs,ref_genome) 
p1 <- plot_spectrum(type_occurrences,by=c("cf","wb"),CT=TRUE)

png("output/figures/type_occ.png",width=5, height=3,units="in",res=500,type="cairo")
p1
dev.off()


###download COSMIC SBS signatures
sbs_url <- "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh38.txt"
dbs_url <- "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3_DBS_GRCh38.txt"
indel_url <- "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3_ID_GRCh37.txt"

sbs <- read.table(sbs_url,header=TRUE)
dbs <- read.table(dbs_url,header=TRUE)
indel <- read.table(indel_url,header=TRUE)

rownames(sbs) <- sbs$Type
rownames(dbs) <- dbs$Type %>% str_replace(., "_", ">")
rownames(indel) <- indel$Type
sbs_sig <- as.matrix(sbs[,2:80])
dbs_sig <- as.matrix(dbs[,2:12])
indel_sig <- as.matrix(indel[,2:19])

cos_sim_samples_signatures = cos_sim_matrix(mut_mat, sbs_sig)
plot_cosine_heatmap(cos_sim_samples_signatures,
                    col_order = sbs_order,
                    cluster_rows = TRUE)

fit_sbs <- fit_to_signatures(mut_mat, sbs_sig)

fit_dbs <- fit_to_signatures(dbs_counts, dbs_sig)
fit_indel <- fit_to_signatures()

# Select signatures with some contribution
select <- which(rowSums(fit_sbs$contribution) > 10)

select_dbs <- which(rowSums(fit_dbs$contribution) > 0.5)
# Plot contribution barplot
plot_contribution(fit_sbs$contribution[select,],
                  sbs_sig[,select],
                  coord_flip = FALSE,
                  mode = "absolute")

plot_contribution(fit_dbs$contribution[select_dbs,],
                  dbs_sig[,select_dbs],
                  coord_flip = FALSE,
                  mode = "relative")

##Indel context 
##create chain (hg38 -> hg19) for liftover
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
ch

##lift over gr.test from hg38 to hg19
gr.indel_19 <- liftOver(grl.comp.indel,ch)
unlist(gr.test_19) -> gr.test_19
genome(gr.indel_19) <- "hg19"


##indel context
plot_indel_contexts(count_indel_contexts(get_indel_context(gr.indel_19,hg19)))
plot_main_indel_contexts(count_indel_contexts(get_indel_context(gr.indel_19,hg19)))
