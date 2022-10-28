# Barchart and Overlap Barplot_2A.jpeg  Barplot_2B.jpeg  Barplot_S1.jpeg  Barplot_S2.jpeg

##Barchart for number of significant genes with p<0.20 FDR cutoff
setwd("/home/workstation/Documents/Sarthak/gender_prognosis")
#for intfile
library(plyr)
library(readr)
library(stringr)
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/")
outfile <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/max_min_reg.txt"
myfiles = list.files(pattern="*_FDR_int.tab", full.names=FALSE)
int_df <- data.frame()
for (i in myfiles){
  int <- read.table(file = i,  sep = '\t', header = TRUE)
  if(nrow(int)>1){
    print(i)
    cancer = sapply(strsplit(i,"_"),"[",1)
    int <- int[int$Interaction_coef == max(int$Interaction_coef)| int$Interaction_coef == min(int$Interaction_coef) ,]
    int_df <- rbind(int_df, data.frame(cancer, int$gene))
  }
}
write.table(int_df,outfile,row.names = F,col.names = F,sep=" ")
