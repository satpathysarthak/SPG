setwd("/home/workstation/Documents/Sarthak/gender_prognosis")
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210122_power")
library(data.table)
library(ggplot2)
library(stringr)
input_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/gender"
input_files  <- list.files(input_dir,pattern = "_int.tab",full.names = TRUE)
old_df <- data.frame()
cut = 0.05
for(i in 1:length(input_files)){
  int <- fread(input_files[i])
  int <- subset(int, select = c(1,24))
  initial <- nrow(int)
  name = sapply(strsplit(input_files[i],split = "/"),'[',9)
  name = sapply(strsplit(name,split = "_"),'[',1)
  int <- int[int$Interaction_Pr...z..< cut,]
  final <- nrow(int)
  new_df <- data.frame("name"=name,"initial"=initial,"final"=final)
  old_df <- rbind(old_df,new_df)
}
data_cutoff <- old_df
rm(list=setdiff(ls(),'data_cutoff'))
data_PFI <- fread("/home/workstation/Documents/Sarthak/gender_prognosis/analysis/table/20210128_list.txt")
data <- merge(data_cutoff,data_PFI, by.x="name",by.y = "cancer")
data_sub <- data[,c(1,5,6)]
library(reshape)
data_sub <- melt(data_sub, id="name")
ggp <- ggplot(data_sub, aes(x = reorder(name,-value), y = value, fill = variable, label = value)) +  # Create stacked bar chart
  geom_bar(stat = "identity")+ xlab("cancer")+ylab("frequency of patients")
ggp
ggp2 <- ggplot(data_cutoff, aes(x = reorder(name,-final), y = final)) +  
  geom_bar(stat = "identity")+ xlab("cancer")+ylab("frequency of genes")
ggp2
ggp3 <- ggplot(data, aes(x = reorder(name,-PFI_1), y = PFI_1)) +  # Create stacked bar chart
  geom_bar(stat = "identity")+ xlab("cancer")+ylab("frequency of PFI_1 patients")
ggp3
data$PFI_ratio = data$PFI_1/data$PFI_0
ggp4 <- ggplot(data, aes(x = reorder(name,-PFI_ratio), y = PFI_ratio)) +  
  geom_bar(stat = "identity")+ xlab("cancer")+ylab("ratio PFI_1:PFI_0")
ggp4
data$PFI_ratio = data$PFI_1/data$patients
ggp5 <- ggplot(data, aes(x = reorder(name,-PFI_ratio), y = PFI_ratio)) +  
  geom_bar(stat = "identity")+ xlab("cancer")+ylab("ratio PFI_1:total")
ggp5
ggp6 <- ggplot(data, aes(x = patients, y = final)) +  
  geom_point(aes(color = factor(name)))+ xlab("patients")+ylab("genes")
ggp6

cor.test(data$final,data$patients,method="spearman")
cor.test(data$final,data$PFI_1,method="spearman")
data$PFI_ratio = data$PFI_1/data$PFI_0
cor.test(data$final,data$PFI_ratio,method="spearman")
data$PFI_ratio = data$PFI_1/data$patients
cor.test(data$final,data$PFI_ratio,method="spearman")
