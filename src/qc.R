# ______________________________________________________________________________
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: QC of samples
#
# Input: seqdata, QC fastq files
#
# Output: QC Sequencing runs, Heatmap, Comparison table between germline profiles
#         of all samples one patient has
# press ALT-O
# ______________________________________________________________________________
#####   Dependencies   #####
library(base)
library(dplyr)
library(xlsx)
library(readxl)

#####   QC sequencing runs  ####
HSM_1346 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "HSMetrics-P1346", 
                      col_types = c("text", "numeric", "numeric"))
GS_1346 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "General Statistics-P1346", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))
HSM_1519 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "HSMetrics-P1519", 
                       col_types = c("text", "numeric", "numeric"))
GS_1519 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "General Statistics-P1519", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))
HSM_1803 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "HSMetrics-P1803", 
                       col_types = c("text", "numeric", "numeric"))
GS_1803 <- read_excel("data/raw/QC/QC_runs.xlsx", sheet = "General Statistics-P1803", 
                      col_types = c("text", "numeric", "numeric", 
                                    "numeric"))


left_join(HSM_1346,GS_1346)->QC_1346
left_join(HSM_1519,GS_1519)->QC_1519
left_join(HSM_1803,GS_1803)->QC_1803


for (file in c("QC_1803","QC_1519","QC_1346")){
  print(file)
  eval(as.name(file))->QC
  #remove samples that didnt work-> "OvCA_45_cf_C1D1.realigned" "OvCA_46_cf_C1D1.realigned" "OvCA_54_cf_C1D1.realigned"
  #QC%>%
    #filter(`Target Bases 30X` <= 0.95)->failed
  #failed$`Sample Name`->failed
  #QC%>%
    #filter(!is.element(`Sample Name`, failed))->QC
  summary(QC)->d
  rm(QC)
  
  filename="output/qc/QC-sequencing runs_means.xlsx"
  write.xlsx(d,filename,sheetName=file,append=TRUE)
  rm(d)
}


#####   Germline Comparison across all samples for each patient   ####
load('data/interim/seqdata.RData')
df %>% 
  filter(snp == TRUE)%>%
  filter(avsnp150!=".")%>%
  filter(FisherScore>1)%>%
  filter(AF >= 0.3)-> df.snp
##Preparation for Heatmap
df.snp%>%
  mutate(PatSample = paste(Patient.ID,Sample,Visite,sep="_"))%>%
  mutate(Sample.ID.match = paste (Sample.ID,replicate,sep="_"))->df.snp



#####      Heatmap ####
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
      
      germline_comparison <- rbind(germline_comparison, data.frame(patient_id = patient_id, sample_id_1 = sample_ids[i], sample_id_2 = sample_ids[j], overlap_percentage = overlap_percentage))
    }
  }
}

# View the results
germline_comparison

filename="output/qc/Distanz.xlsx"
write.xlsx(germline_comparison,filename,sheetName="germline comparison",append=TRUE)

rm(list=ls())


