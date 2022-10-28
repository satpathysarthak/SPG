input_dir <- "/media/sarthak/gender_prognosis/cox_regression/gender"
input_files  <- list.files(input_dir,pattern = "_int.tab",full.names = TRUE)
setwd("/media/sarthak/gender_prognosis/plots/20201112_pval_hist")
library(data.table)
library(ggplot2)
library(stringr)
for(i in 1:length(input_files)){
  int <- fread(input_files[i])
  cancer <- strsplit(input_files[i], split = "/", fixed = T)[[1]]
  cancer = cancer[length(cancer)]
  cancer <- strsplit(cancer, split = "_", fixed = T)[[1]]
  cancer = cancer[1]
  #int$cancer = cancer
  #int <- int[int$Interaction_Pr...z..<0.5,]
  #int$p_adjusted_new <- p.adjust(int$Interaction_Pr...z.., method = "BH")
  #int <- int[int$p_adjusted_new < 0.10,]
  # ggplot(aes(log10(Interaction_Pr...z..)), data = int)+
  #   geom_histogram(binwidth = 30)+
  #   ggtitle(cancer)
  #dev.off()
  jpeg(filename = str_c(cancer,"_int.jpeg"))
  hist(log10(int$Interaction_Pr...z..),xlab = "log(P_int)",main=cancer)
  dev.off()
  print(cancer)
}
