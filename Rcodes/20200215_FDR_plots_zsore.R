## FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.20
library(data.table)
library(stringr)
setwd("/home/workstation/Documents/Sarthak/gender_prognosis")
##### Gender ####
#print("Gender")
input_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/gender_zscore"
output_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore"
input_files  <- list.files(input_dir,pattern = "_est.tab",full.names = TRUE)
output_files <- sub(input_dir,output_dir,input_files)
output_files <- sub("_est.tab","_FDR_est.tab",output_files)
#print(input_dir)
old_df <- data.frame()
for(i in 1:length(input_files)){
  est <- fread(input_files[i])
  a<- nrow(est)
  #print(str_c(input_files[i],nrow(est),ncol(est),sep = "\t"))
  est <- est[abs(est$gene_coef) > quantile(abs(est$gene_coef),0.8),]
  #shapiro.test(est$gene_coef)
  #print(nrow(est))
  b<- nrow(est)
  est$p_adjusted_new <- p.adjust(est$gene_Pr...z.., method = "BH")
  est <- est[est$p_adjusted_new < 0.10,]
  c<- nrow(est)
  ##write.table(est,output_files[i],append = FALSE, col.names = TRUE,row.names = FALSE, sep='\t',quote=FALSE)
  #print(str_c(output_files[i],nrow(est),ncol(est),sep = "\t"))
  new_df <- data.frame("file"=input_files[i], 'initial'= a,"intermediate"=b,"final"=c)
  old_df <- rbind(old_df,new_df)
  rm(list="est")
  #print("##")
}
rm(list=c("input_files","output_files","i"))
#print("est to int")
input_files  <- list.files(input_dir,pattern = "_int.tab",full.names = TRUE)
output_files <- sub(input_dir,output_dir,input_files)
output_files <- sub("_int.tab","_FDR_int.tab",output_files)
for(i in 1:length(input_files)){
  int <- fread(input_files[i])
  a<- nrow(int)
  #print(str_c(input_files[i],nrow(int),ncol(int),sep = "\t"))
  hist(abs(int$Interaction_coef))
  print(shapiro.test(sample(abs(int$Interaction_coef),size=500)))
  int <- int[abs(int$Interaction_coef) > quantile(abs(int$Interaction_coef),0.8),]
  b<- nrow(int)
  #print(nrow(int))
  int$p_adjusted_new <- p.adjust(int$Interaction_Pr...z.., method = "BH")
  int <- int[int$p_adjusted_new < 0.10,]
  c<- nrow(int)
  #write.table(int,output_files[i],append = FALSE, col.names = TRUE,row.names = FALSE, sep='\t',quote=FALSE)
  #print(str_c(output_files[i],nrow(int),ncol(int),sep = "\t"))
  new_df <- data.frame("file"=input_files[i],"initial"= a,"intermediate"=b,"final"=c)
  old_df <- rbind(old_df,new_df)
  rm(list="int")
  #print("##")
}
#rm(list = ls())

######### miR Gender
setwd("/home/workstation/Documents/Sarthak/gender_prognosis")
#print("##############################")
#print("miR Gender")
input_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/miR_gender"
output_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_miR_gender"
input_files  <- list.files(input_dir,pattern = "_est.tab",full.names = TRUE)
output_files <- sub(input_dir,output_dir,input_files)
output_files <- sub("_est.tab","_FDR_est.tab",output_files)
#print(input_dir)
for(i in 1:length(input_files)){
  est <- fread(input_files[i])
  a<- nrow(est)
  #print(str_c(input_files[i],nrow(est),ncol(est),sep = "\t"))
  est <- est[abs(est$gene_coef)>mean(abs(est$gene_coef)),]
  #print(nrow(est))
  b<- nrow(est)
  est$p_adjusted_new <- p.adjust(est$gene_Pr...z.., method = "BH")
  est <- est[est$p_adjusted_new < 0.20,]
  c<- nrow(est)
  ##write.table(est,output_files[i],append = FALSE, col.names = TRUE,row.names = FALSE, sep='\t',quote=FALSE)
  #print(str_c(output_files[i],nrow(est),ncol(est),sep = "\t"))
  new_df <- data.frame("file"=input_files[i],"initial"= a,"intermediate"=b,"final"=c)
  old_df <- rbind(old_df,new_df)
  rm(list="est")
  #print("##")
}
rm(list=c("input_files","output_files","i"))
#print("est to int")
input_files  <- list.files(input_dir,pattern = "_int.tab",full.names = TRUE)
output_files <- sub(input_dir,output_dir,input_files)
output_files <- sub("_int.tab","_FDR_int.tab",output_files)
for(i in 1:length(input_files)){
  int <- fread(input_files[i])
  a<- nrow(int)
  #print(str_c(input_files[i],nrow(int),ncol(int),sep = "\t"))
  int <- int[abs(int$gene_coef)>mean(abs(int$gene_coef)),]
  b<- nrow(int)
  #print(nrow(int))
  int$p_adjusted_new <- p.adjust(int$gene_Pr...z.., method = "BH")
  int <- int[int$p_adjusted_new < 0.20,]
  c<- nrow(int)
  #write.table(int,output_files[i],append = FALSE, col.names = TRUE,row.names = FALSE, sep='\t',quote=FALSE)
  #print(str_c(output_files[i],nrow(int),ncol(int),sep = "\t"))
  new_df <- data.frame("file"=input_files[i],"initial"= a,"intermediate"=b,"final"=c)
  old_df <- rbind(old_df,new_df)
  rm(list="int")
  #print("##")
}
data <- old_df
colnames(data) <- c("file_name","median_1fpkm","mean_absolute","BH")
data$file_name <- as.character(data$file_name)
data$file_name <- sub("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/","",data$file_name)
data$type <- sapply(strsplit(data$file_name, split = "/",fixed = TRUE),"[",1)
data$cancer <- sapply(strsplit(data$file_name, split = "/",fixed = TRUE),"[",2)
data$model <- sapply(strsplit(data$cancer, split = "_",fixed = TRUE),"[",2)
data$model <- sub(".tab","",data$model)
data$cancer <- sapply(strsplit(data$cancer, split = "_",fixed = TRUE),"[",1)
data <- subset(data, select = c("type","model","cancer","median_1fpkm","mean_absolute","BH"))
colnames(data)[5] <- "quantile"
#write.table(data,"/home/workstation/Documents/Sarthak/gender_prognosis/analysis/table/20200813_FDR_filters.txt",quote = FALSE, sep = "\t",row.names=FALSE)
#### Plotting ####

setwd("/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210214_FDR_exploratory_zscore")
library(reshape2)
library('tidyverse')
library(ggplot2)
results_long <- melt(data, id.vars=c("cancer","model","type"))
jpeg(filename = "gender_est.jpeg")
results_long1 <- results_long[results_long$model == "est"&results_long$type=="gender_zscore",]
results_long1 %>% 
  #filter(type %in% variable) %>%
  ggplot(aes(x=cancer,y=value, fill=variable))+
  geom_bar(stat="identity", position=position_dodge()) + xlab("cancer type") + ylab("frequency")+
  ggtitle("gender_est")
dev.off()
jpeg(filename = "gender_int.jpeg")
results_long1 <- results_long[results_long$model == "int"&results_long$type=="gender_zscore",]
results_long1 %>% 
  #filter(type %in% variable) %>%
  ggplot(aes(x=cancer,y=value, fill=variable))+
  geom_bar(stat="identity", position=position_dodge()) + xlab("cancer type") + ylab("frequency")+
  ggtitle("gender_int")
dev.off()
jpeg(filename = "miR_gender_est.jpeg")
results_long1 <- results_long[results_long$model == "est"&results_long$type=="miR_gender",]
results_long1 %>% 
  #filter(type %in% variable) %>%
  ggplot(aes(x=cancer,y=value, fill=variable))+
  geom_bar(stat="identity", position=position_dodge()) + xlab("cancer type") + ylab("frequency")+
  ggtitle("miR_gender_est")
dev.off()
jpeg(filename = "miR_gender_int.jpeg")
results_long1 <- results_long[results_long$model == "int"&results_long$type=="miR_gender",]
results_long1 %>% 
  #filter(type %in% variable) %>%
  ggplot(aes(x=cancer,y=value, fill=variable))+
  geom_bar(stat="identity", position=position_dodge()) + xlab("cancer type") + ylab("frequency")+
  ggtitle("miR_gender_int")
dev.off()


#est vs int
jpeg(filename = "FDR_gender.jpeg")
results_long1 <- results_long[results_long$variable == "BH"&results_long$type=="gender_zscore",]
results_long1 %>% 
  #filter(type %in% variable) %>%
  ggplot(aes(x=cancer,y=value, fill=model))+
  geom_bar(stat="identity", position=position_dodge()) + xlab("cancer type") + ylab("frequency")+
  ggtitle("FDR")
dev.off()
