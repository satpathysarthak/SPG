# Genes in the category: ns_ns,s_ns,ns_s,s_s
library(data.table)
library(ggplot2)
input_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/bivariate_gender/"
input_files_male  <- list.files(input_dir,pattern = "_male.tab",full.names = F)
input_files_female  <- list.files(input_dir,pattern = "_female.tab",full.names = F)
map_dir = "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/withnames/"
input_files_map <- list.files(map_dir,pattern = "_map_int.txt",full.names = T)
outdir <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/"
setwd(input_dir)
data_old <- data.frame()
for(i in 1:length(input_files_male)){
  map <- fread(input_files_map[i], header = F)
  colnames(map) <- c("gene","gene_name")
  data_male <- fread(input_files_male[i])
  data_male <- data_male[,c(1:6)]
  data_female <- fread(input_files_female[i])
  data_female <- data_female[,c(1:6)]
  colnames(data_male)[c(2:6)] <- gsub("gene","male_gene",colnames(data_male)[c(2:6)])
  colnames(data_female)[c(2:6)] <- gsub("gene","female_gene",colnames(data_female)[c(2:6)])
  data <- merge(data_male,data_female,by='gene')
  data$direction_male <- ifelse(data$male_gene_Pr...z..<0.05,1,0)
  data$direction_female <- ifelse(data$female_gene_Pr...z..<0.05,1,0)
  print(table(data$direction_male,data$direction_female))
  cancer <- sapply(strsplit(input_files_male[i], split = "_",fixed = TRUE),"[",1)
  dir <- as.vector(table(data$direction_male,data$direction_female))
  data$sig <- NA
  data$sig <- as.character(paste0(data$direction_male,data$direction_female))
  data$sig_dir <- data$male_gene_coef*data$female_gene_coef
  data$sig_dir <- data$sig_dir/abs(data$sig_dir)
  data <- merge(data,map, by='gene')
  ns_ns <- nrow(data[data$sig == "00",])
  s_ns <- nrow(data[data$sig == "10",])
  ns_s <- nrow(data[data$sig == "01",])
  s_s <- nrow(data[data$sig == "11",])
  s_s_opp <- subset(data, sig == "11" & sig_dir == "-1")
  s_s_same <- subset(data, sig == "11" & sig_dir == "1")
  data_old <- rbind(data_old,data.frame(cancer,ns_ns,s_ns,ns_s,s_s))
  write.table(data[data$sig == "00",], paste0(outdir,cancer,"_ns_ns",".tab"), append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
  write.table(data[data$sig == "10",], paste0(outdir,cancer,"_s_ns",".tab"), append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
  write.table(data[data$sig == "01",], paste0(outdir,cancer,"_ns_s",".tab"), append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
  #write.table(data[data$sig == "11",], paste0(outdir,cancer,"_s_s",".tab"), append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
  write.table(s_s_opp, paste0(outdir,cancer,"_s_s_opp",".tab"), append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
  write.table(s_s_same, paste0(outdir,cancer,"_s_s_same",".tab"), append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
}
