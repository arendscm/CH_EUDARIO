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
library(base)
library(dplyr)
library(xlsx)
library(stringr)
library(reshape)
library(tidyr)
library(readxl)
library(reshape2)


########   set working directory #####
#setwd('H:/Meine Ablage')
#setwd("C:/Users/maxar/Documents/AG Damm/EUDARIO/data_analysis/EUDARIO")

########  Load preprocessed sequencing data
#df <- read.csv('data/interim/mutationcalls.csv')
load('data/interim/seqdata.RData')

######## Get Patient ids
source("src/ids.R")



##SNP
df %>% 
  filter(snp == TRUE)%>%
  filter(avsnp150!=".")%>%
  filter(FisherScore>1)%>%
  filter(AF >= 0.3)-> df.snp
##Preparation for Heatmap
df.snp%>%
  mutate(PatSample = paste(Patient.ID,Sample,Visite,sep="_"))%>%
  mutate(Sample.ID.match = paste (Sample.ID,replicate,sep="_"))->df.snp

########   SNP-Comparison -> Heatmap ####
#QC_germline<-read.table(file = 'QC/Total_result.txt', sep = '\t', header = FALSE)
#names(QC_germline) <- c("Sample1","Sample2","Score", "MatchingStatus")
#filename="QC/QC_germline.xlsx"
#write.xlsx(QC_germline,filename, sheetName="QC", append=TRUE)


##rescue matching germline Mutations
semi_join(df, df.snp, by="Patmut")->df.snp

#Heatmap (Robert Skript)
### long table to wide table
newdata <- dcast(data = df.snp.Pat,
                 formula = position ~ Sample.ID.match,    
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
pdf(file="output/qc/Distanznewdata.pdf",height=35,width=45)
image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, rownames(dst), cex.axis = 0.5, las=3)
axis(2, 1:dim, rownames(dst), cex.axis = 0.5, las=1)

text(expand.grid(1:dim, 1:dim), sprintf("%1i", dst) , cex=0.6)
dev.off()

#write newdata in excelfile
filename="output/qc/Distanz.xlsx"
write.xlsx(newdata,filename,sheetName="Distanz",append=TRUE)

rm(newdata)
rm(dim)
rm(dst)
rm(filename)

##comparison germline:
# Get the unique patient IDs
patient_ids <- unique(df.snp$Patient.ID)
sample_ids <- unique(df.snp$Sample.ID)



# Initialize a data frame to store the results
germline_comparison <- data.frame(patient_id = character(0), sample_id_1 = character(0), sample_id_2 = character(0), overlap_percentage = numeric(0))

# Loop through each patient ID
for (patient_id in patient_ids) {
  
  # Subset the data frame to only include samples for the current patient
  patient_samples <- df.snp[df.snp$Patient.ID == patient_id, ]
  
  # Get the unique sample IDs for the patient
  sample_ids <- unique(patient_samples$Sample.ID.match)
  
  # If there is only one sample for the patient, skip the patient
  if (length(sample_ids) < 2) {
    next
  }
  
  # Loop through each sample combination
  for (i in 1:(length(sample_ids)-1)) {
    for (j in (i+1):length(sample_ids)) {
      
      # Subset the data frame to only include the current samples
      sample_1 <- patient_samples[patient_samples$Sample.ID.match== sample_ids[i], ]
      sample_2 <- patient_samples[patient_samples$Sample.ID.match == sample_ids[j], ]
      
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

filename="output/qc/Distanz.xlsx"
write.xlsx(germline_comparison,filename,sheetName="germline comparison",append=TRUE)

rm(list=ls())

