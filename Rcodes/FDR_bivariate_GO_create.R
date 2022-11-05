# Creates a ranked file from output of FunrichR
output_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/"
dir_path <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/"
library(stringr)
setwd(dir_path)
dir <- list.files(dir_path, pattern= ".tab", all.files=FALSE,  full.names=FALSE)
for(file_name in dir){
  print(file_name)
  id_file <- file_name
  id <- read.delim(id_file, header = TRUE)
  id <- subset(id, select = c('gene_name'))
  if(nrow(id)<1){next}
  output_file <- strsplit(file_name, split = ".",fixed=T)[[1]][1]
  output_file <- str_c(output_dir,output_file,'.rnk')
  id$gene_name <- as.character(id$gene_name)
  id <- unique.data.frame(id)
  #write.table('# rnk format file', file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)
  write.table(id, file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append = F)
}
  
