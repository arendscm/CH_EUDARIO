library(base)
library(dplyr)
library(stringr)
library(reshape)
library(tidyr)
library(reshape2)


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

########   Create txt.file for Circleplot ####--------------------------------------------

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
write.table(finaltable, file = "output/Circleplot.txt", sep = " ",
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
