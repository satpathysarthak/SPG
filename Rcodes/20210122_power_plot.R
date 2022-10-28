setwd("/media/sarthak/gender_prognosis")
setwd("/media/sarthak/gender_prognosis/plots/20210122_power")
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
  #int$order <- log10(int$Interaction_Pr...z..)
  #cut = c(ceiling(min(int$order)):0)
  cut = c(0.00005,0.0001,0.0005, 0.001, 0.005,0.01,0.05,0.1,0.2)
  power = data.frame("cutoff"=cut)
  power$number = NA
  for(j in 1:length(cut)){
    power$number[j] = nrow(int[int$Interaction_Pr...z..<cut[j]])
  }
  #png(filename = str_c(name,".png"))
  plot = ggplot(data = power,aes(x=log10(cutoff),y=number))+
    geom_line(stat="identity") + xlab("log10(cutoff)") + ylab("frequency")+
    geom_bar(stat="identity", position=position_dodge()) + xlab("log10(cutoff)") + ylab("frequency")+
    ggtitle(str_c(name,' ',nrow(int),' genes'))
  #geom_hline(yintercept=0.05*abs(min(cut))/nrow(int),color = "red")+
  #geom_hline(yintercept=0.05*abs(min(cut))/nrow(int[int$Interaction_Pr...z..<0.05]),color = "blue")
  ggsave(plot, file=paste0("20200122_", name,".png"), width = 44.45, height = 27.78, units = "cm", dpi=300)
    #dev.off()
}
