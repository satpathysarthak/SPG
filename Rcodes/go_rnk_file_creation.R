#ID mapping: Specific from taking input from .int file
####
input_file_ref_from <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_ensembl.tab"
input_file_ref_to <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_Gene_Name.tab"
#id_file <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis_FDR/intfile2/GBM_int.tab2"
output_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/go/"
dir_path <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/"
###
library(stringr)
setwd(dir_path)
dir <- list.files(dir_path, pattern= "_FDR_int.tab", all.files=FALSE,  full.names=FALSE)
dir <- dir[-c(2,10,11,12)]
#dir <- "THCA_int.tab2"
for(file_name in dir){
  print(file_name)
  id_file <- file_name
  ref <- read.delim(input_file_ref_from, header=FALSE)
  id <- read.table(id_file, quote="\"", comment.char="", header = TRUE)
  id <- subset(id, select = c('gene', "Interaction_Pr...z.."))
  ref$V1 <- as.character(ref$V1)
  ref$V3 <- as.character(ref$V3)
  id$gene <- as.character(id$gene)
  id$gene <- sapply(strsplit(as.character(id$gene), "\\."), `[`, 1)
  for(i in 1:nrow(id)){
    for(j in 1:nrow(ref)){
      if(id$gene[i] == ref$V3[j]){
        id$gene_uni[i] = ref$V1[j]
        break
      }
    }
  }
  ref <- read.delim(input_file_ref_to, header=FALSE)
  ref$V1 <- as.character(ref$V1)
  ref$V3 <- as.character(ref$V3)
  id$gene_uni <- as.character(id$gene_uni)
  #id$V3 <- as.character(id$V3)
  for(i in 1:nrow(id)){
    for(j in 1:nrow(ref)){
      if(id$gene_uni[i] == ref$V1[j]){
        id$gene_name[i] = ref$V3[j]
        break
        
      }
    }
  }
  id <- subset(id, select = c('gene_name', 'Interaction_Pr...z..'))
  id <- id[order(-id$Interaction_Pr...z..),]
  output_file <- strsplit(file_name, split = "_")[[1]][1]
  output_file <- str_c(output_dir,output_file,'.rnk')
  write.table('# rnk format file', file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)
  write.table(id, file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
}
