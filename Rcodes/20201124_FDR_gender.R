## FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.20
library(data.table)
library(stringr)
setwd("/media/sarthak/gender_prognosis")
##### Gender ####
print("Gender")
input_dir <- "/media/sarthak/gender_prognosis/cox_regression/gender"
output_dir <- "/media/sarthak/gender_prognosis/cox_regression/20201124_1_FDR_gender"
input_files  <- list.files(input_dir,pattern = "_est.tab",full.names = TRUE)
output_files <- sub(input_dir,output_dir,input_files)
output_files <- sub("_est.tab","_FDR_est.tab",output_files)
print(input_dir)
for(i in 1:length(input_files)){
  est <- fread(input_files[i])
  print(str_c(input_files[i],nrow(est),ncol(est),sep = "\t"))
  est <- est[est$gene_Pr...z..<0.05,]
  print(nrow(est))
  est$p_adjusted_new <- p.adjust(est$gene_Pr...z.., method = "BH")
  est <- est[est$p_adjusted_new < 0.10,]
  write.table(est,output_files[i],append = FALSE, col.names = TRUE,row.names = FALSE, sep='\t',quote=FALSE)
  print(str_c(output_files[i],nrow(est),ncol(est),sep = "\t"))
  rm(list="est")
  print("##")
}
rm(list=c("input_files","output_files","i"))
print("est to int")
input_files  <- list.files(input_dir,pattern = "_int.tab",full.names = TRUE)
output_files <- sub(input_dir,output_dir,input_files)
output_files <- sub("_int.tab","_FDR_int.tab",output_files)
for(i in 1:length(input_files)){
  int <- fread(input_files[i])
  print(str_c(input_files[i],nrow(int),ncol(int),sep = "\t"))
  int <- int[int$Interaction_Pr...z..<0.05,]
  print(nrow(int))
  int$p_adjusted_new <- p.adjust(int$Interaction_Pr...z.., method = "BH")
  int <- int[int$p_adjusted_new < 0.05,]
  write.table(int,output_files[i],append = FALSE, col.names = TRUE,row.names = FALSE, sep='\t',quote=FALSE)
  print(str_c(output_files[i],nrow(int),ncol(int),sep = "\t"))
  rm(list="int")
  print("##")
}
rm(list = ls())
