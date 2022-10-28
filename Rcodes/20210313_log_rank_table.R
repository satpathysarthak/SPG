library(data.table)

input_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/log_rank/"
input_files  <- list.files(input_dir,pattern = ".name",full.names = F)
setwd(input_dir)
data_old <- data.frame()
for(i in 1:length(input_files)){
  log_data <- fread(input_files[i])
  log_data$direction_male <- ifelse(log_data$male_p<0.05,1,0)
  log_data$direction_female <- ifelse(log_data$female_p<0.05,1,0)
  print(table(log_data$direction_male,log_data$direction_female))
  cancer <- sapply(strsplit(input_files[i], split = "_",fixed = TRUE),"[",1)
  dir <- as.vector(table(log_data$direction_male,log_data$direction_female))
  ns_ns <- dir[1]
  s_ns <- dir[2]
  ns_s <- dir[3]
  s_s <- dir[4]
  data_old <- rbind(data_old,data.frame(cancer,ns_ns,s_ns,ns_s,s_s))
}
