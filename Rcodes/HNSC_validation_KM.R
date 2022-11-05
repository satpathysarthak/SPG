# KM  plots for HNSC validation and training data

#TCGA data

library(data.table)
library(stringr)
library(dplyr)
library(survival)
library(survivalROC)

args <- 'HNSC'
print(args)
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/")

risk_calc <- function(dataframe,coefficient_df){
  risk = vector(length=nrow(dataframe))
  for(i in 1:length(coefficient_df)){
    risk = risk + (dataframe[names(coefficient_df)[i]]*coefficient_df[names(coefficient_df)[i]])
  }
  names(risk)= NULL
  #row.names(risk) = NULL
  return(risk)
}

#### Input & Output file names ####
#tumor_file = str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Tumor/",args[1],".tab")
clinical_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Clinical/processed/",args[1],"_pro.tab")
normalized_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/data/",args[1],".tab")
model_boot_both_file <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/HNSC_SPG_PG_both.tab' 
model_boot_est_file <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/HNSC_SPG_PG_est.tab' 
#output_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/",args[1],"_SPG.tab")
#####initialising Sample Dataset#####
sample <- fread(clinical_file)
sample <- unique.data.frame(sample)
sample$ID <- sapply(strsplit(sample$submitter_id,split = "-"), "[",3)
sample$gender <- as.character(sample$gender)
sample[sample$gender=="female",]$gender <- "1"
sample[sample$gender=="male",]$gender <- "0"
sample$gender <- as.numeric(sample$gender)
sample <- sample[!is.na(sample$gender),]
##### Data Merging #####
data_org <- read.table(normalized_file,sep='\t',header = T, stringsAsFactors = F)
model_both <- fread(model_boot_both_file)
model_est <- fread(model_boot_est_file)
model_no <- 91
cox_prog_both <- as.numeric(strsplit(model_both$coefficients[model_both$bootstrap==model_no],';')[[1]])
names(cox_prog_both) <- strsplit(model_both$predictors[model_both$bootstrap==model_no],';')[[1]]
cox_prog_est <- as.numeric(strsplit(model_est$coefficients[model_est$bootstrap==model_no],';')[[1]])
names(cox_prog_est) <- strsplit(model_est$predictors[model_est$bootstrap==model_no],';')[[1]]


data = data_org
data = data_org[data_org$gender == 0,]
data = data_org[data_org$gender == 1,]
data$Risk_Score_both <- risk_calc(data,cox_prog_both)
data$Risk_Score_est <- risk_calc(data,cox_prog_est)
data$Risk_Score_both <- unlist(data$Risk_Score_both)
data$Risk_Score_est <- unlist(data$Risk_Score_est)
med <- median(data$Risk_Score_both)
# med = 0.00075947540994856 # all
# med = -0.0629065300017545 # male
# med = 0.135008114710481 # female
#HNSC
# med = 0.008299481 # all
# med = 0.04970602 # male
# med = -0.03123197 # female

data['Risk_Score_both']<-ifelse(data['Risk_Score_both']< med, 0, 1)



rm(list = ls())

## GEO data
input_file <- '/home/workstation/Documents/Sarthak/gender_prognosis/validation/HNSC/output/data_HNSC.tab'
model_boot_both_file <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/HNSC_SPG_PG_both.tab' 
model_boot_est_file <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/HNSC_SPG_PG_est.tab' 
map_genes <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/HNSC_model_genes_mapped_91.txt'

data <- read.table(file = input_file,  sep = '\t', header = TRUE, stringsAsFactors = F)
data_org = data
model_both <- fread(model_boot_both_file)
model_est <- fread(model_boot_est_file)
map <- fread(map_genes)


risk_calc <- function(dataframe,coefficient_df){
  coefficient_df <- coefficient_df[names(coefficient_df)%in%names(dataframe)]
  risk = vector(length=nrow(dataframe))
  for(i in 1:length(coefficient_df)){
    risk = risk + (dataframe[names(coefficient_df)[i]]*coefficient_df[names(coefficient_df)[i]])
  }
  names(risk)= NULL
  #row.names(risk) = NULL
  return(risk)
}

##
model_no <- 91
cox_prog_both <- as.numeric(strsplit(model_both$coefficients[model_both$bootstrap==model_no],';')[[1]])
names(cox_prog_both) <- strsplit(model_both$predictors[model_both$bootstrap==model_no],';')[[1]]
cox_prog_est <- as.numeric(strsplit(model_est$coefficients[model_est$bootstrap==model_no],';')[[1]])
names(cox_prog_est) <- strsplit(model_est$predictors[model_est$bootstrap==model_no],';')[[1]]

df_both_genes <- data.frame('gene' = names(cox_prog_both),'coef'=cox_prog_both)
df_est_genes <- data.frame('gene' = names(cox_prog_est), 'coef'=cox_prog_est)

df_both_genes <- merge(df_both_genes, map, by='gene', all.y=F)
df_est_genes <- merge(df_est_genes, map, by='gene', all.y=F)

df_both_genes <- unique.data.frame(df_both_genes)
df_est_genes <- unique.data.frame(df_est_genes)

df_both_genes <- df_both_genes[df_both_genes$gene_name != '',]
df_est_genes <- df_est_genes[df_est_genes$gene_name != '',]

cox_prog_both <- df_both_genes$coef
names(cox_prog_both) <- df_both_genes$gene_name
cox_prog_est <- df_est_genes$coef
names(cox_prog_est) <- df_est_genes$gene_name

#data = data_org[data_org$gender==1,]
data = data_org
data = data_org[data_org$gender==0,]
data = data_org[data_org$gender==1,]
data$Risk_Score_both <- risk_calc(data,cox_prog_both)
data$Risk_Score_est <- risk_calc(data,cox_prog_est)
data$Risk_Score_both <- unlist(data$Risk_Score_both)
data$Risk_Score_est <- unlist(data$Risk_Score_est)
med <- median(data$Risk_Score_both)
# med = 0.008299481 # all
# med = 0.04970602 # male
# med = -0.03123197 # female
data['Risk_Score_both']<-ifelse(data['Risk_Score_both']< med, 0, 1)
data$Risk_Score_both <- unlist(data$Risk_Score_both)
fit <- survfit( Surv(time = pfs,
                     event = pfs_event)~ unlist(data['Risk_Score_both']),
                data =data)
ggsurvplot(fit, data=data,palette = c("#404040", "#ca0020"),
           ggtheme = theme_bw(),pval = TRUE,
           pval.method=TRUE,conf.int = FALSE,
           font.main=c(20,"bold"),font.x=c(16,"bold"),font.y=c(16,"bold"),
           font.tickslab=c(14,"plain"), title = 'GEO combined model: ', legend.labs = c("Low", "High"))

data = data_org
foo <- data[data$gender==1,]
foo2 <- data[data$gender==0,]
foo2 <- foo2[sample (c(1:nrow(foo2)), size=47, replace =F),]
data <- rbind(foo,foo2)
