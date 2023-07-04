# ______________________________________________________________________________
# Ovarian Cancer filtering Script
#
# Author: Max & Klara
#
# Description: global definition of interesting gene groups
#
# Input: none
#
# Output: gene group definitions
#
# ______________________________________________________________________________
#####  Dependencies   #####

##interesting gene groups
ch_genes <- c("DNMT3A","TET2","ASXL1","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5","U2AF1","TP53","PPM1D","RAD21","PTPN11","KIT","BCORL1","BCOR","ATM")
tp53_genes <- c("TP53")
ppm1d_genes <- c("PPM1D")
ddr_genes <- c("TP53","PPM1D","ATM","CHEK2","ATR")
brca_genes <- c("BRCA1","BRCA2")
hrd_genes <- c("BRCA1","BRCA2","ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")
ch_genes_without_HRD <- c("DNMT3A","TET2","ASXL1","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5","U2AF1","TP53","PPM1D","RAD21","PTPN11","KIT","BCORL1","BCOR")
ovca_genes <- c("TP53","NF1","CDK12","BRCA1","BRCA2")



