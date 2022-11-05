# Creates a scatter plot male and female coefficients 
library(data.table)
library(ggplot2)
input_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/bivariate_gender/"
input_files_male  <- list.files(input_dir,pattern = "_male.tab",full.names = F)
input_files_female  <- list.files(input_dir,pattern = "_female.tab",full.names = F)
outdir <- "/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210413_male_female_int_plots/"
int_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/withnames"
input_int <- list.files(int_dir,pattern = "_int.tab",full.names = T)
setwd(input_dir)
data_old <- data.frame()
for(i in 1:length(input_files_male)){
  data_int <- fread(input_int[i])
  data_male <- fread(input_files_male[i])
  data_male <- data_male[,c(1:6)]
  data_female <- fread(input_files_female[i])
  data_female <- data_female[,c(1:6)]
  colnames(data_male)[c(2:6)] <- gsub("gene","male_gene",colnames(data_male)[c(2:6)])
  colnames(data_female)[c(2:6)] <- gsub("gene","female_gene",colnames(data_female)[c(2:6)])
  data <- merge(data_male,data_female, by='gene')
  data  <-merge(data,data_int, by='gene')
  data$direction_male <- ifelse(data$male_gene_Pr...z..<0.05,"s","ns")
  data$direction_female <- ifelse(data$female_gene_Pr...z..<0.05,"s","ns")
  print(table(data$direction_male,data$direction_female))
  cancer <- sapply(strsplit(input_files_male[i], split = "_",fixed = TRUE),"[",1)
  per_male <- quantile(abs(data$male_gene_coef),0.50)
  per_female <- quantile(abs(data$female_gene_coef),0.50)
  data$Interaction_log_p <- -log10(data$Interaction_Pr...z..)
  
  colorset = c('ns'='#404040','s'='#6fbde3')
  sp1 <- ggplot(data = data, aes(x=Interaction_coef,y=male_gene_coef,color=direction_male))+
    geom_point()+ xlab(paste0("Interaction Beta")) + ylab("Male Beta") + ggtitle(paste0(cancer,"_male"))
  # sp <- sp+scale_color_manual(values=c("#404040","#F1B6DA","#92C5DE","red", "darkgrey"))
  sp1 <- sp1+scale_color_manual(values=colorset)
  sp1 <- sp1+theme_classic()+geom_hline(yintercept=0,linetype = "dotted")+geom_vline(xintercept=0,linetype = "dotted")
  #sp <- sp + xlim(c(lab_all*-1,lab_all)) + ylim(c(lab_all*-1,lab_all))
  sp1 <- sp1 + guides(color = guide_legend(override.aes = list(size = 0.5)))
  #colorset = c('ns_ns'='#404040','ns_s'='#F1B6DA','s_ns'='#92C5DE','s_s_opp'='green','s_s_same'='darkgrey')
  sp1 <- sp1 + geom_text_repel(aes(label=ifelse((abs(male_gene_coef)>per_male) ,as.character(gene_name),'')),size=2)
  tiff(paste0(outdir,cancer,"_int_beta_male.tiff"),height=1200,width=1834,res=200)
  print(sp1)
  dev.off()
  
  colorset = c('ns'='#404040','s'='#f77ec8')
  sp2 <- ggplot(data = data, aes(x=Interaction_coef,y=female_gene_coef,color=direction_female))+
    geom_point()+ xlab(paste0("Interaction Beta")) + ylab("Female Beta") + ggtitle(paste0(cancer,"_female"))
  # sp <- sp+scale_color_manual(values=c("#404040","#F1B6DA","#92C5DE","red", "darkgrey"))
  sp2 <- sp2+scale_color_manual(values=colorset)
  sp2 <- sp2+theme_classic()+geom_hline(yintercept=0,linetype = "dotted")+geom_vline(xintercept=0,linetype = "dotted")
  #sp <- sp + xlim(c(lab_all*-1,lab_all)) + ylim(c(lab_all*-1,lab_all))
  sp2 <- sp2 + guides(color = guide_legend(override.aes = list(size = 0.5)))
  #colorset = c('ns_ns'='#404040','ns_s'='#F1B6DA','s_ns'='#92C5DE','s_s_opp'='green','s_s_same'='darkgrey')
  sp2 <- sp2 + geom_text_repel(aes(label=ifelse((abs(female_gene_coef)>per_female) ,as.character(gene_name),'')),size=2)
  tiff(paste0(outdir,cancer,"_int_beta_female.tiff"),height=1200,width=1834,res=200)
  print(sp2)
  dev.off()
}

rm(list=ls())

#### In this bock size is by the p-val
library(data.table)
library(ggplot2)
input_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/bivariate_gender/"
input_files_male  <- list.files(input_dir,pattern = "_male.tab",full.names = F)
input_files_female  <- list.files(input_dir,pattern = "_female.tab",full.names = F)
outdir <- "/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210413_male_female_int_plots/"
int_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/withnames"
input_int <- list.files(int_dir,pattern = "_int.tab",full.names = T)
setwd(input_dir)
data_old <- data.frame()
for(i in 1:length(input_files_male)){
  data_int <- fread(input_int[i])
  data_male <- fread(input_files_male[i])
  data_male <- data_male[,c(1:6)]
  data_female <- fread(input_files_female[i])
  data_female <- data_female[,c(1:6)]
  colnames(data_male)[c(2:6)] <- gsub("gene","male_gene",colnames(data_male)[c(2:6)])
  colnames(data_female)[c(2:6)] <- gsub("gene","female_gene",colnames(data_female)[c(2:6)])
  data <- merge(data_male,data_female, by='gene')
  data  <-merge(data,data_int, by='gene')
  data$direction_male <- ifelse(data$male_gene_Pr...z..<0.05,"s","ns")
  data$direction_female <- ifelse(data$female_gene_Pr...z..<0.05,"s","ns")
  print(table(data$direction_male,data$direction_female))
  cancer <- sapply(strsplit(input_files_male[i], split = "_",fixed = TRUE),"[",1)
  per_male <- quantile(abs(data$male_gene_coef),0.50)
  per_female <- quantile(abs(data$female_gene_coef),0.50)
  data$Interaction_log_p <- -log10(data$Interaction_Pr...z..)
  
  colorset = c('ns'='#404040','s'='#6fbde3')
  sp1 <- ggplot(data = data, aes(x=Interaction_coef,y=male_gene_coef,color=direction_male))+
    geom_point(aes(size=Interaction_log_p))+ xlab(paste0("Interaction Beta")) + ylab("Male Beta") + ggtitle(paste0(cancer,"_male"))
  # sp <- sp+scale_color_manual(values=c("#404040","#F1B6DA","#92C5DE","red", "darkgrey"))
  sp1 <- sp1 + scale_size(range = c(0, 10)) + labs(size="-log10(Int p-value)",col="Male Beta p-val")
  sp1 <- sp1+scale_color_manual(values=colorset)
  sp1 <- sp1+theme_classic()+geom_hline(yintercept=0,linetype = "dotted")+geom_vline(xintercept=0,linetype = "dotted")
  #sp <- sp + xlim(c(lab_all*-1,lab_all)) + ylim(c(lab_all*-1,lab_all))
  sp1 <- sp1 + guides(color = guide_legend(override.aes = list(size = 0.5)))
  #colorset = c('ns_ns'='#404040','ns_s'='#F1B6DA','s_ns'='#92C5DE','s_s_opp'='green','s_s_same'='darkgrey')
  sp1 <- sp1 + geom_text_repel(aes(label=ifelse((abs(male_gene_coef)>per_male) ,as.character(gene_name),'')),size=2)
  tiff(paste0(outdir,cancer,"_int_beta_male.tiff"),height=1200,width=1834,res=200)
  print(sp1)
  dev.off()
  
  colorset = c('ns'='#404040','s'='#f77ec8')
  sp2 <- ggplot(data = data, aes(x=Interaction_coef,y=female_gene_coef,color=direction_female))+
    geom_point(aes(size=Interaction_log_p))+ xlab(paste0("Interaction Beta")) + ylab("Female Beta") + ggtitle(paste0(cancer,"_female"))
  # sp <- sp+scale_color_manual(values=c("#404040","#F1B6DA","#92C5DE","red", "darkgrey"))
  sp2 <- sp2 + scale_size(range = c(0, 10)) + labs(size="-log10(Int p-value)",col="Female Beta p-val")
  sp2 <- sp2+scale_color_manual(values=colorset)
  sp2 <- sp2+theme_classic()+geom_hline(yintercept=0,linetype = "dotted")+geom_vline(xintercept=0,linetype = "dotted")
  #sp <- sp + xlim(c(lab_all*-1,lab_all)) + ylim(c(lab_all*-1,lab_all))
  sp2 <- sp2 + guides(color = guide_legend(override.aes = list(size = 0.5)))
  #colorset = c('ns_ns'='#404040','ns_s'='#F1B6DA','s_ns'='#92C5DE','s_s_opp'='green','s_s_same'='darkgrey')
  sp2 <- sp2 + geom_text_repel(aes(label=ifelse((abs(female_gene_coef)>per_female) ,as.character(gene_name),'')),size=2)
  tiff(paste0(outdir,cancer,"_int_beta_female.tiff"),height=1200,width=1834,res=200)
  print(sp2)
  dev.off()
}  


library(data.table)
library(ggplot2)
input_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/bivariate_gender/"
input_files_male  <- list.files(input_dir,pattern = "_male.tab",full.names = F)
input_files_female  <- list.files(input_dir,pattern = "_female.tab",full.names = F)
outdir <- "/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210413_male_female_int_plots/"
int_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/withnames"
input_int <- list.files(int_dir,pattern = "_int.tab",full.names = T)
setwd(input_dir)
data_old <- data.frame()
for(i in 1:length(input_files_male)){
  data_int <- fread(input_int[i])
  data_male <- fread(input_files_male[i])
  data_male <- data_male[,c(1:6)]
  data_female <- fread(input_files_female[i])
  data_female <- data_female[,c(1:6)]
  colnames(data_male)[c(2:6)] <- gsub("gene","male_gene",colnames(data_male)[c(2:6)])
  colnames(data_female)[c(2:6)] <- gsub("gene","female_gene",colnames(data_female)[c(2:6)])
  data <- merge(data_male,data_female, by='gene')
  data  <-merge(data,data_int, by='gene')
  data$direction_male <- ifelse(data$male_gene_Pr...z..<0.05,"s","ns")
  data$direction_female <- ifelse(data$female_gene_Pr...z..<0.05,"s","ns")
  print(table(data$direction_male,data$direction_female))
  cancer <- sapply(strsplit(input_files_male[i], split = "_",fixed = TRUE),"[",1)
  per_male <- quantile(abs(data$male_gene_coef),0.50)
  per_female <- quantile(abs(data$female_gene_coef),0.50)
  data$Interaction_log_p <- -log10(data$Interaction_Pr...z..)
  
  colorset = c('ns'='#404040','s'='#6fbde3')
  sp1 <- ggplot(data = data, aes(x=Interaction_coef,y=male_gene_coef,color=direction_male))+
    geom_point()+ xlab(paste0("Interaction Beta")) + ylab("Male Beta") + ggtitle(paste0(cancer,"_male"))
  # sp <- sp+scale_color_manual(values=c("#404040","#F1B6DA","#92C5DE","red", "darkgrey"))
  sp1 <- sp1+scale_color_manual(values=colorset)
  sp1 <- sp1+theme_classic()+geom_hline(yintercept=0,linetype = "dotted")+geom_vline(xintercept=0,linetype = "dotted")
  #sp <- sp + xlim(c(lab_all*-1,lab_all)) + ylim(c(lab_all*-1,lab_all))
  sp1 <- sp1 + guides(color = guide_legend(override.aes = list(size = 0.5)))
  #colorset = c('ns_ns'='#404040','ns_s'='#F1B6DA','s_ns'='#92C5DE','s_s_opp'='green','s_s_same'='darkgrey')
  sp1 <- sp1 + geom_text_repel(aes(label=ifelse((abs(male_gene_coef)>per_male) ,as.character(gene_name),'')),size=2)
  tiff(paste0(outdir,cancer,"_int_beta_male.tiff"),height=1200,width=1834,res=200)
  print(sp1)
  dev.off()
  
  colorset = c('ns'='#404040','s'='#f77ec8')
  sp2 <- ggplot(data = data, aes(x=Interaction_coef,y=female_gene_coef,color=direction_female))+
    geom_point()+ xlab(paste0("Interaction Beta")) + ylab("Female Beta") + ggtitle(paste0(cancer,"_female"))
  # sp <- sp+scale_color_manual(values=c("#404040","#F1B6DA","#92C5DE","red", "darkgrey"))
  sp2 <- sp2+scale_color_manual(values=colorset)
  sp2 <- sp2+theme_classic()+geom_hline(yintercept=0,linetype = "dotted")+geom_vline(xintercept=0,linetype = "dotted")
  #sp <- sp + xlim(c(lab_all*-1,lab_all)) + ylim(c(lab_all*-1,lab_all))
  sp2 <- sp2 + guides(color = guide_legend(override.aes = list(size = 0.5)))
  #colorset = c('ns_ns'='#404040','ns_s'='#F1B6DA','s_ns'='#92C5DE','s_s_opp'='green','s_s_same'='darkgrey')
  sp2 <- sp2 + geom_text_repel(aes(label=ifelse((abs(female_gene_coef)>per_female) ,as.character(gene_name),'')),size=2)
  tiff(paste0(outdir,cancer,"_int_beta_female.tiff"),height=1200,width=1834,res=200)
  print(sp2)
  dev.off()
}

rm(list=ls())

library(data.table)
library(ggplot2)
library(viridis)
input_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/bivariate_gender/"
input_files_male  <- list.files(input_dir,pattern = "_male.tab",full.names = F)
input_files_female  <- list.files(input_dir,pattern = "_female.tab",full.names = F)
outdir <- "/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210413_male_female_int_plots/"
int_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/withnames"
input_int <- list.files(int_dir,pattern = "_int.tab",full.names = T)
setwd(input_dir)
data_old <- data.frame()
for(i in 1:length(input_files_male)){
  data_int <- fread(input_int[i])
  data_male <- fread(input_files_male[i])
  data_male <- data_male[,c(1:6)]
  data_female <- fread(input_files_female[i])
  data_female <- data_female[,c(1:6)]
  colnames(data_male)[c(2:6)] <- gsub("gene","male_gene",colnames(data_male)[c(2:6)])
  colnames(data_female)[c(2:6)] <- gsub("gene","female_gene",colnames(data_female)[c(2:6)])
  data <- merge(data_male,data_female, by='gene')
  data  <-merge(data,data_int, by='gene')
  data$direction_male <- ifelse(data$male_gene_Pr...z..<0.05,"s","ns")
  data$direction_female <- ifelse(data$female_gene_Pr...z..<0.05,"s","ns")
  print(table(data$direction_male,data$direction_female))
  cancer <- sapply(strsplit(input_files_male[i], split = "_",fixed = TRUE),"[",1)
  per_male <- quantile(abs(data$male_gene_coef),0.50)
  per_female <- quantile(abs(data$female_gene_coef),0.50)
  data$Interaction_log_p <- -log10(data$Interaction_Pr...z..)
  
  shapeset = c('ns'=2,'s'=1)
  sp1 <- ggplot(data = data, aes(x=Interaction_coef,y=male_gene_coef,shape=direction_male, color=Interaction_log_p))+
    geom_point()+ xlab(paste0("Interaction Beta")) + ylab("Male Beta") + ggtitle(paste0(cancer,"_male"))
  # sp <- sp+scale_color_manual(values=c("#404040","#F1B6DA","#92C5DE","red", "darkgrey"))
  sp1 <- sp1+scale_shape_manual(values=shapeset)
  sp1 <- sp1+scale_color_viridis(option = "B", direction = -1)
  sp1 <- sp1+theme_classic()+geom_hline(yintercept=0,linetype = "dotted")+geom_vline(xintercept=0,linetype = "dotted")
  #sp <- sp + xlim(c(lab_all*-1,lab_all)) + ylim(c(lab_all*-1,lab_all))
  sp1 <- sp1 + guides(color = guide_legend(override.aes = list(size = 0.5)))
  #colorset = c('ns_ns'='#404040','ns_s'='#F1B6DA','s_ns'='#92C5DE','s_s_opp'='green','s_s_same'='darkgrey')
  sp1 <- sp1 + geom_text_repel(aes(label=ifelse((abs(male_gene_coef)>per_male) ,as.character(gene_name),'')),size=2)
  tiff(paste0(outdir,cancer,"_int_beta_male.tiff"),height=1200,width=1834,res=200)
  print(sp1)
  dev.off()
  
  shapeset = c('ns'=2,'s'=1)
  sp2 <- ggplot(data = data, aes(x=Interaction_coef,y=female_gene_coef,shape=direction_female, color=Interaction_log_p))+
    geom_point()+ xlab(paste0("Interaction Beta")) + ylab("Female Beta") + ggtitle(paste0(cancer,"_female"))
  # sp <- sp+scale_color_manual(values=c("#404040","#F1B6DA","#92C5DE","red", "darkgrey"))
  sp2 <- sp2+scale_shape_manual(values=shapeset)
  sp2 <- sp2+scale_color_viridis(option = "B", direction = -1)
  sp2 <- sp2+theme_classic()+geom_hline(yintercept=0,linetype = "dotted")+geom_vline(xintercept=0,linetype = "dotted")
  #sp <- sp + xlim(c(lab_all*-1,lab_all)) + ylim(c(lab_all*-1,lab_all))
  sp2 <- sp2 + guides(color = guide_legend(override.aes = list(size = 0.5)))
  #colorset = c('ns_ns'='#404040','ns_s'='#F1B6DA','s_ns'='#92C5DE','s_s_opp'='green','s_s_same'='darkgrey')
  sp2 <- sp2 + geom_text_repel(aes(label=ifelse((abs(female_gene_coef)>per_female) ,as.character(gene_name),'')),size=2)
  tiff(paste0(outdir,cancer,"_int_beta_female.tiff"),height=1200,width=1834,res=200)
  print(sp2)
  dev.off()
}


