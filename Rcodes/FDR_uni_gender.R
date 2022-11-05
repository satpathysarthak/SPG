## FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.20
library(data.table)
library(stringr)
setwd("/media/sarthak/gender_prognosis")
##### Gender ####
print("Gender")
input_dir <- "/media/sarthak/gender_prognosis/cox_regression/univariate_gender"
output_dir <- "/media/sarthak/gender_prognosis/cox_regression/FDR_univariate_gender"
input_files  <- list.files(input_dir,pattern = "_male.tab",full.names = TRUE)
output_files <- sub(input_dir,output_dir,input_files)
output_files <- sub("_male.tab","_FDR_male.tab",output_files)
print(input_dir)
for(i in 1:length(input_files)){
  male <- fread(input_files[i])
  print(str_c(input_files[i],nrow(male),ncol(male),sep = "\t"))
  male <- male[abs(male$coef)>mean(abs(male$coef)),]
  print(nrow(male))
  male$p_adjusted_new <- p.adjust(male$Pr...z.., method = "BH")
  male <- male[male$p_adjusted_new < 0.20,]
  write.table(male,output_files[i],append = FALSE, col.names = TRUE,row.names = FALSE, sep='\t',quote=FALSE)
  print(str_c(output_files[i],nrow(male),ncol(male),sep = "\t"))
  rm(list="male")
  print("##")
}
rm(list=c("input_files","output_files","i"))
print("male to female")
input_files  <- list.files(input_dir,pattern = "_female.tab",full.names = TRUE)
output_files <- sub(input_dir,output_dir,input_files)
output_files <- sub("_female.tab","_FDR_female.tab",output_files)
for(i in 1:length(input_files)){
  female <- fread(input_files[i])
  print(str_c(input_files[i],nrow(female),ncol(female),sep = "\t"))
  female <- female[abs(female$coef)>mean(abs(female$coef)),]
  print(nrow(female))
  female$p_adjusted_new <- p.adjust(female$Pr...z.., method = "BH")
  female <- female[female$p_adjusted_new < 0.20,]
  write.table(female,output_files[i],append = FALSE, col.names = TRUE,row.names = FALSE, sep='\t',quote=FALSE)
  print(str_c(output_files[i],nrow(female),ncol(female),sep = "\t"))
  rm(list="female")
  print("##")
}
rm(list = ls())
