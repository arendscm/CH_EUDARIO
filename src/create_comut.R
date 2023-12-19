# ______________________________________________________________________________
# CH in EUDARIO
#
# Author: Max & Klara
#
# Description: definition of createComut function, along the lines of 
# Bernard et al. Nature Communications volume 12, Article number: 7244 (2021)
#
# Input: seqdata
#
# Output: comutation plot
#
# ______________________________________________________________________________
#####  Dependencies   #####
library(base)
library(dplyr)
library(forcats)

createComut <- function(sub,dis){

# calculate frequency of gene 1 occuring before gene 2
genes <- as.data.frame(unique(sub$Gene))
results_list <- list()
n=1

for(i in 1:nrow(genes)){
  print(i)
  # i<- 1
  # select mutations of interest
  gene_1 <- as.character(genes[i,])
  i2=(i+1)
  if(i2>nrow(genes)){next}
  for(j in i2:nrow(genes)){
    #print(j)
    # j <- 2
    gene_2 <- as.character(genes[j,])
    
    if(gene_1 != gene_2){
      sub_gene_1 <- subset(sub, sub$Gene == gene_1)
      sub_gene_2 <- subset(sub, sub$Gene == gene_2)
      
      # select patients with both mutations
      gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
      
      n_cases <- as.numeric(nrow(gene_1_and_2))
      
      # add point color columns for visualizing clonal/subclonal trends
      gene_1_and_2$vaf_ratio <- as.numeric((gene_1_and_2$TVAF.x - gene_1_and_2$TVAF.y))
      print(as.numeric((gene_1_and_2$TVAF.x - gene_1_and_2$TVAF.y)))
      # gene_1_and_2$vaf_ratio <- as.numeric(gene_1_and_2$vaf_ratio)
      
      gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
      
      if(nrow(gene_1_and_2) > 0){
        # define point color
        gene_1_and_2$Clonality <- NA
        
        for(k in 1:nrow(gene_1_and_2)){
          if(gene_1_and_2$vaf_ratio[k] <= dis & gene_1_and_2$vaf_ratio[k] >= -dis){
            gene_1_and_2$Clonality[k] <- 0
          }
          if(gene_1_and_2$vaf_ratio[k] > dis){
            gene_1_and_2$Clonality[k] <- 1
          }
          if(gene_1_and_2$vaf_ratio[k] < -dis){
            gene_1_and_2$Clonality[k] <- 2
          }
        }
        
        gene_1_and_2$Clonality <- as.numeric(gene_1_and_2$Clonality)
        
        # fraction of cases where gene x occurs before gene y
        n_1_before_2 <- as.numeric(length(which(gene_1_and_2$Clonality == 1)))
        n_2_before_1 <- as.numeric(length(which(gene_1_and_2$Clonality == 2)))
        
        temp_dat <- data.frame(matrix(NA, nrow = 1, ncol = 4))
        colnames(temp_dat) = c("Gene_1", "Gene_2", "number_1_before_2", "number_2_before_1")
        
        temp_dat[n,1] <- (gene_1) 
        temp_dat[n,2] <- (gene_2)
        temp_dat[n,3] <- (n_1_before_2) 
        temp_dat[n,4] <- (n_2_before_1) 
        
        results_list[[n]] <- temp_dat
        n=n+1     
      }
    }
  }
}
temp_final = na.omit(as.data.frame(do.call(rbind, results_list)))

temp = temp_final[,c(2,1,3,4)]
names(temp) = c("Gene_1", "Gene_2", "number_2_before_1", "number_1_before_2")

temp_final = rbind(temp_final, temp)


temp_final$fraction_1_then_2 = ((temp_final$number_1_before_2)/((temp_final$number_1_before_2)+(temp_final$number_2_before_1))*100)


# now find the number of co-occuring cases
genes <- unique(sub$Gene)
temp_dat_2 <- data.frame(matrix(NA, nrow = length(genes), ncol = length(genes)))

rownames(temp_dat_2) <- genes
colnames(temp_dat_2) <- genes

for(i in 1:nrow(temp_dat_2)){
  # select mutations of interest
  gene_x <- rownames(temp_dat_2)[i]
  
  for(j in 1:ncol(temp_dat_2)){
    gene_y <- colnames(temp_dat_2)[j]
    
    if(gene_x != gene_y){
      sub_gene_1 <- subset(sub, sub$Gene == gene_x)
      sub_gene_2 <- subset(sub, sub$Gene == gene_y)
      
      # select patients with both mutations
      gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
      
      if(nrow(gene_1_and_2) > 0){
        # add point color columns for visualizing clonal/subclonal trends
        gene_1_and_2$vaf_ratio <- (gene_1_and_2$TVAF.x - gene_1_and_2$TVAF.y)
        
        gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$Sample), ] 
        
        # number of cases where gene x and gene y co-occur
        temp_dat_2[i,j] <- as.numeric(nrow(gene_1_and_2))
      }  
    }
  }
}


# Get lower triangle of the correlation matrix
get_upper_tri<-function(temp_dat_2){
  temp_dat_2[lower.tri(temp_dat_2)] <- NA
  return(temp_dat_2)
}

temp_dat_2_final <- get_upper_tri(temp_dat_2)

temp_dat_2_final$genes <- rownames(temp_dat_2_final)

temp_dat_2_final_melted_1 <- melt(temp_dat_2_final, na.rm = TRUE)

colnames(temp_dat_2_final_melted_1) <- c("Gene_1", "Gene_2", "number_1_and_2")

# Get lower triangle of the correlation matrix
get_lower_tri<-function(temp_dat_2){
  temp_dat_2[upper.tri(temp_dat_2)] <- NA
  return(temp_dat_2)
}

temp_dat_2_final <- get_lower_tri(temp_dat_2)

temp_dat_2_final$genes <- rownames(temp_dat_2_final)

temp_dat_2_final_melted_2 <- melt(temp_dat_2_final, na.rm = TRUE)

colnames(temp_dat_2_final_melted_2) <- c("Gene_1", "Gene_2", "number_1_and_2")

temp_dat_final_melted_2 <- rbind(temp_dat_2_final_melted_1, temp_dat_2_final_melted_2)

#combine number and fraction of supporting cases
temp_dat_final_melted <- left_join(temp_dat_final_melted_2, temp_final, by = c("Gene_1", "Gene_2"))

temp_dat_2_final_melted <- temp_dat_final_melted[,c(2,1,3)]
colnames(temp_dat_2_final_melted) <- c("Gene_1", "Gene_2", "number_1_and_2")

temp_dat_final_melted <- left_join(temp_dat_final_melted, temp_dat_2_final_melted, by = c("Gene_1", "Gene_2","number_1_and_2"))

# temp_dat_final_melted <- temp_dat_final_melted[order(temp_dat_final_melted$Gene_2),]


temp_dat_final_melted$bin <- NA

for(i in 1:nrow(temp_dat_final_melted)){
  if(temp_dat_final_melted$number_1_and_2[i] <= 2){
    temp_dat_final_melted$bin[i] <- "2"
  }
  if(temp_dat_final_melted$number_1_and_2[i] > 2 & temp_dat_final_melted$number_1_and_2[i] <= 4){
    temp_dat_final_melted$bin[i] <- "2-4"
  }
  if(temp_dat_final_melted$number_1_and_2[i] > 4 & temp_dat_final_melted$number_1_and_2[i] <= 9){
    temp_dat_final_melted$bin[i] <- "5-9"
  }
  if(temp_dat_final_melted$number_1_and_2[i] > 9 & temp_dat_final_melted$number_1_and_2[i] <= 20){
    temp_dat_final_melted$bin[i] <- "10-20"
  }
  if(temp_dat_final_melted$number_1_and_2[i] > 20){
    temp_dat_final_melted$bin[i] <- ">20"
  }
}

temp_dat_final_melted$Gene_2 <- as.character(temp_dat_final_melted$Gene_2)
return(temp_dat_final_melted)
}

