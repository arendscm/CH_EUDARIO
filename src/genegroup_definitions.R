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
typical_ch_genes <- c("DNMT3A","TET2","ASXL1","CBL","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5","U2AF1","TP53","PPM1D","RAD21","PTPN11","ATM","STAG2","CHEK2","RUNX1","STAT3","BCOR","BCORL1","BRCC3")
ch_genes <- c("DNMT3A","TET2","ASXL1","CBL","CEBPA","GNB1","GNAS","IDH1","IDH2","JAK2","SF3B1","SRSF2","U2AF1;U2AF1L5","U2AF1","TP53","PPM1D","CHEK2","RAD21","PTPN11","KIT","BCORL1","BCOR","ATM","STAG2","STAT3","KRAS","NRAS","ETV6","CSF3R","WT1","SETBP1","MYD88","FLT3","NF1","RUNX1","EZH2","GATA1","GATA2","CALR","MPL","NPM1","PHF6","BRCC3","BRAF","NOTCH1","XPO1")
tp53_genes <- c("TP53")
ppm1d_genes <- c("PPM1D")
ddr_genes <- c("TP53","PPM1D","ATM","CHEK2","ATR")
brca_genes <- c("BRCA1","BRCA2")
hrd_genes <- c("BRCA1","BRCA2","ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","EMSY","FAM175A","FANCA","FANCC","FANCI","FANCL","MLH1","MRE11","MSH2","MSH6","NBN","PALB2","PMS2","RAD21","RAD50","RAD51","RAD51C","RAD51D","RAD52","RAD54L","PTEN","BRCC3")
ch_genes_without_HRD <- setdiff(ch_genes,hrd_genes)
typical_ch_genes_without_HRD <- setdiff(typical_ch_genes,hrd_genes)
ovca_genes <- c("TP53","NF1","CDK12","BRCA1","BRCA2")



