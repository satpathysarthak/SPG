setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore")

outdir <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/comp/est_int_overlap/"
myfiles1 = list.files(pattern="*_FDR_int.tab", full.names=FALSE)
myfiles2 = list.files(pattern="*_FDR_est.tab", full.names=FALSE)
venn_df <- data.frame()
for(i in 1:length(myfiles1)){
  int <- read.table(file = myfiles1[i],  sep = '\t', header = TRUE)
  int <- subset(int, select = 'gene')
  int$gene <- as.character(int$gene)
  est <- read.table(file = myfiles2[i],  sep = '\t', header = TRUE)
  est <- subset(est, select = 'gene')
  est$gene <- as.character(est$gene)
  Cancer <- strsplit(myfiles1[i],split="_")[[1]][1]
  print(Cancer)
  print(strsplit(myfiles2[i],split="_")[[1]][1])
  Int <- length(int$gene)
  Est <- length(est$gene)
  Overlap <- length(intersect(est$gene, int$gene))
  v <- cbind(Cancer,Int,Est,Overlap)
  venn_df = rbind(venn_df,v)
  write.table(intersect(est$gene, int$gene), str_c(outdir,Cancer,"_overlap.txt"), col.names = FALSE
              , row.names = FALSE, quote = FALSE, sep = "\t")
  print("#####")
}
