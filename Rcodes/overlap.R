## Calculating overlap with X, Y chromosome; gtex (tissue specific); sex-biased genes 


# chromosome X Y
setwd("/home/workstation/Documents/Sarthak/gender_prognosis")
#for intfile
library(plyr)
library(readr)
library(stringr)
library(data.table)
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/")
myfiles = list.files(pattern="*_FDR_int.tab", full.names=FALSE)
chrX_file = '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/chr/chr_X_genes.tsv'
chrY_file = '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/chr/chr_Y_genes.tsv'
chrX = fread(chrX_file, header = F, stringsAsFactors = F)
chrY = fread(chrY_file, header = F,stringsAsFactors = F)
int_df <- data.frame()
for (i in myfiles){
  int <- read.table(file = i,  sep = '\t', header = TRUE,stringsAsFactors = F)
  tot <- nrow(int)
  if(tot <1){next ;}
  int$ens = sapply(strsplit(int$gene, ".", fixed = TRUE), "[",1)
  print <- strsplit(i,split="_")[[1]][1]
  noX <- length(intersect(chrX$V3, int$ens))
  noY <- length(intersect(chrY$V3, int$ens))
  int_df <- rbind(int_df, data.frame("Cancer"=print, "Int"=tot, 'X'=noX, 'Y'=noY))
}

#for estfile
myfiles = list.files(pattern="*_FDR_est.tab", full.names=FALSE)
est_df <- data.frame()
for (i in myfiles){
  est <- read.table(file = i,  sep = '\t', header = TRUE, stringsAsFactors = F)
  tot <- nrow(est)
  if(tot <1){next ;}
  est$ens = sapply(strsplit(est$gene, ".", fixed = TRUE), "[",1)
  print <- strsplit(i,split="_")[[1]][1]
  noX <- length(intersect(chrX$V3, est$ens))
  noY <- length(intersect(chrY$V3, est$ens))
  est_df <- rbind(est_df, data.frame("Cancer"=print, "Est"=tot, 'X'=noX, 'Y'=noY))
}
rm(list = setdiff(ls(),c("est_df","int_df")))
outfile_int = '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/chr/overlap_int.txt'
outfile_est = '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/chr/overlap_est.txt'
write.table(int_df, outfile_int, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
write.table(est_df, outfile_est, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
rm(list=ls())

####
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/withnames/")
myfiles = list.files(pattern="*_FDR_name_int.tab", full.names=FALSE)
cell_dir = '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/cell_pap/mRNA_miRNA_sex_bias'
cell_files = list.files(path = cell_dir,pattern="*_mRNA_sex_bias.txt", full.names=T)
int_df <- data.frame()
for (i in myfiles){
  int <- read.table(file = i,  sep = '\t', header = TRUE,stringsAsFactors = F)
  tot <- nrow(int)
  if(tot <1){next ;}
  print <- strsplit(i,split="_")[[1]][1]
  cell <- fread(grep(print, cell_files, value = T),
                sep = '|', header = F, stringsAsFactors = F)
  #int$ens = sapply(strsplit(int$gene, ".", fixed = TRUE), "[",1)
  noCell <- length(intersect(cell$V1, int$gene_name))
  #noY <- length(intersect(chrY$V3, int$ens))
  int_df <- rbind(int_df, data.frame("Cancer"=print, "Int"=tot, 'Cell'=noCell))
}
outfile_int = '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/cell_pap/overlap_int.txt'
write.table(int_df, outfile_int, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
rm(list=ls())

## 

gtex_file = '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/gtex/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt'
gtex <- fread(gtex_file, header = T, stringsAsFactors = F)
gtex$gene <- sapply(strsplit(gtex$gene, ".", fixed = TRUE), "[",1)
gtex <- subset(gtex, select = c('gene','tissue'))
#gtex <-  gtex[, `:=` (count = .N), by = gene]
#gtex <- subset(gtex, select = c('gene','gene'))
gtex <- unique.data.frame(gtex)

setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/")
myfiles = list.files(pattern="*_FDR_int.tab", full.names=FALSE)
map_file <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/gtex/map_tissue.txt'
map <- fread(map_file, stringsAsFactors = F)
int_df <- data.frame()
for (i in myfiles){
  int <- read.table(file = i,  sep = '\t', header = TRUE,stringsAsFactors = F)
  tot <- nrow(int)
  if(tot <1){next ;}
  int$ens = sapply(strsplit(int$gene, ".", fixed = TRUE), "[",1)
  print <- strsplit(i,split="_")[[1]][1]
  list_tis <- unique(grep(map[map$cancer == print]$tissue,gtex$tissue, value = T))
  sub_gtex <- gtex[gtex$tissue %in% list_tis]
  sub_gtex <- sub_gtex[, `:=` (count = .N), by = gene]
  sub_gtex <- subset(sub_gtex, select = c('gene','count'))
  sub_gtex <- unique.data.frame(sub_gtex)
  noGtex <- length(intersect(unique(sub_gtex$gene), int$ens))
  print(str_c(i,list_tis, collapse=';'))
  int_df <- rbind(int_df, data.frame("Cancer"=print, "Int"=tot, 'gtex'=noGtex))
}
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/")
myfiles = list.files(pattern="*_FDR_est.tab", full.names=FALSE)
est_df <- data.frame()
for (i in myfiles){
  est <- read.table(file = i,  sep = '\t', header = TRUE,stringsAsFactors = F)
  tot <- nrow(est)
  if(tot <1){next ;}
  est$ens = sapply(strsplit(est$gene, ".", fixed = TRUE), "[",1)
  print <- strsplit(i,split="_")[[1]][1]
  list_tis <- unique(grep(map[map$cancer == print]$tissue,gtex$tissue, value = T))
  sub_gtex <- gtex[gtex$tissue %in% list_tis]
  sub_gtex <- sub_gtex[, `:=` (count = .N), by = gene]
  sub_gtex <- subset(sub_gtex, select = c('gene','count'))
  sub_gtex <- unique.data.frame(sub_gtex)
  noGtex <- length(intersect(unique(sub_gtex$gene), est$ens))
  print(str_c(i,list_tis, collapse=';'))
  est_df <- rbind(est_df, data.frame("Cancer"=print, "Est"=tot, 'gtex'=noGtex))
}

rm(list = setdiff(ls(),c("est_df","int_df")))
outfile_int = '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/gtex/overlap_int.txt'
outfile_est = '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/gtex/overlap_est.txt'
write.table(int_df, outfile_int, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
write.table(est_df, outfile_est, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
rm(list=ls())


