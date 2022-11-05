#ID mapping: Specific from taking input from .int file
####
input_file_ref_from <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_ensembl.tab"
input_file_ref_to <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_Gene_Name.tab"
#id_file <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis_FDR/intfile2/GBM_int.tab2"
#output_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/go/"
dir_path <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/"
###
library(stringr)
setwd(dir_path)
dir <- list.files(dir_path, pattern= "_FDR_int.tab", all.files=FALSE,  full.names=FALSE)
ref <- read.delim(input_file_ref_from, header=FALSE)
id_file = "HNSC_FDR_int.tab"
id <- read.table(id_file, quote="\"", comment.char="", header = TRUE,stringsAsFactors = F)
id <- subset(id, select = c('gene'))
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
genes <- read.delim("/home/workstation/Documents/Sarthak/gender_prognosis/analysis/HNSC_enrich.txt", header=FALSE)
genes <- as.character(unique(genes$V1))
id <- id[id$gene_name %in% genes,]
genes <- id$gene
rm(list = setdiff(ls(),"genes"))
