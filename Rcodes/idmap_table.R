input_file_ref_from <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_ensembl.tab"
input_file_ref_to <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_Gene_Name.tab"
ref <- read.delim(input_file_ref_from, header=FALSE)
id <- read.csv("~/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore_name/a", sep="")
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
