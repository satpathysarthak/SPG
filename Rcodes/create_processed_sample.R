## Merge information on gender with Clinical-OS: Creating processed Clinical Sample File
rm(list = ls())
library(data.table)
library(stringr)
#### listing files and directories
setwd("/media/sarthak/gender_prognosis/")
dir_gender <- "/media/sarthak/gender_prognosis/Clinical/Clinical_cdc/"
files_gender_all <- list.files(dir_gender,pattern = "_gender.tab")
cancer <- sapply(strsplit(files_gender_all, split = "_"),'[',1)
dir_OS <- "/media/sarthak/gender_prognosis/Clinical/Clinical_OS/"
files_OS_all <- list.files(dir_OS)
files_OS <- vector()
for(i in 1:length(cancer)){
  index = grep(cancer[i],files_OS_all)
  files_OS[i] <- files_OS_all[index]
}
rm(list = "files_OS_all","index","i")
grep_list <- c("submitter_id","PFI","PFI.time","^age","^gender$")
dir_out <- "/media/sarthak/gender_prognosis/Clinical/processed/"
####

for(i in 1:length(cancer)){
  print(cancer[i])
  input_file1 <- str_c(dir_OS,files_OS[i])
  print(input_file1)
  file_OS <- fread(input_file1)
  colnames(file_OS)[grep("V1",colnames(file_OS))] <- "submitter_id"
  
#  files_gender <- list.files(dir_gender,pattern = "_gender.tab")
  input_file2 <- str_c(dir_gender,files_gender_all[i])
  files_gender <- fread(input_file2)
  print(input_file2)
  files_gender <- unique.data.frame(files_gender)
  #colnames(file_OS)[grep("submitter_id",colnames(file_OS))] <- "submitter_id"
  file_merged <- merge.data.frame(file_OS,files_gender,by.x = "submitter_id",by.y = "submitter_id", all.x = TRUE, all.y = FALSE)
  print(colnames(file_merged))
  list_colnames <- vector()
  for(l in 1:length(grep_list)){
    list_colnames[l] <- colnames(file_merged)[grep(grep_list[l],colnames(file_merged))]
  }
  print(list_colnames)
  file_merged <- subset(file_merged, select = list_colnames)
  print(colnames(file_merged))
  colnames(file_merged)[grep("^age",colnames(file_merged))] <- "age"
  print(file_merged)
  output_file_name <- str_c(dir_out,cancer[i],"_pro.tab")
  write.table(file_merged,file=output_file_name,quote=F,sep = "\t",row.names = FALSE)
  print(output_file_name)
  rm(list = c("input_file1","input_file2"))
}
warnings()

age_change <- c("KIRP_pro.tab","LUAD_pro.tab","THCA_pro.tab")
for(i in age_change){
  print(i)
  output_file_name <- str_c(dir_out,i)
  file_age <- fread(output_file_name)
  file_age$age <- as.numeric(file_age$age)
  file_age$age <- trunc(file_age$age/365)
  write.table(file_age,file=output_file_name,quote=F,sep = "\t",row.names = FALSE,append = FALSE)
}
