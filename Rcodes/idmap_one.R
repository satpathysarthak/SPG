idmap <- function(name_ensembl){
  name_ensembl <- sapply(strsplit(name_ensembl,".",fixed=T),"[",1)
  input_file_ref_from <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_ensembl.tab"
  input_file_ref_to <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_Gene_Name.tab"
  ref <- read.delim(input_file_ref_from, header=FALSE,stringsAsFactors = F)
  name_uni = "NA"
  for(j in 1:nrow(ref)){
    if(ref$V3[j] == name_ensembl){
      name_uni = ref$V1[j]
      break
    }
  }
  ref <- read.delim(input_file_ref_to, header=FALSE,stringsAsFactors = F)
  name_gene = "NA"
  for(j in 1:nrow(ref)){
    if(ref$V1[j]==name_uni){
      name_gene = ref$V3[j]
      break
      
    }
  }
  if(name_gene == "NA"){name_gene = NA}
  return(name_gene)
}
save(idmap, file = "/home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/idmap_fun.RData")
