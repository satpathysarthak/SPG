setwd("/media/sarthak/gender_prognosis")
setwd("/media/sarthak/gender_prognosis/plots/20210111_power")
library(data.table)
library(ggplot2)
library(stringr)
input_dir <- "/media/sarthak/gender_prognosis/cox_regression/gender"
input_files  <- list.files(input_dir,pattern = "_int.tab",full.names = TRUE)
for(i in 1:length(input_files)){
  int <- fread(input_files[i])
  int <- subset(int, select = c(1,24))
  name = sapply(strsplit(input_files[i],split = "/"),'[',7)
  name = sapply(strsplit(name,split = "_"),'[',1)
  int$order <- log10(int$Interaction_Pr...z..)
  cut = c(ceiling(min(int$order)):0)
  power = data.frame("logP"=cut)
  power$number = NA

  for(j in 1:length(cut)){
    power$number[j] = nrow(int[int$order<cut[j]])
  }
  png(filename = str_c(name,".png"))
  ggplot(data = power,aes(x=logP,y=number))+
    geom_bar(stat="identity", position=position_dodge()) + xlab("logP") + ylab("frequency")+
    ggtitle(name)+
    geom_hline(yintercept=0.05*abs(min(cut))/nrow(int),color = "red")+
    geom_hline(yintercept=0.05*abs(min(cut))/nrow(int[int$Interaction_Pr...z..<0.05]),color = "blue")
  dev.off()
}
