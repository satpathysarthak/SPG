## Script to generate plots for Gene Expression
rm(list = ls())
args<-commandArgs(TRUE)
print(args)
setwd("/media/sarthak/gender_prognosis/")
library(data.table)
library(stringr)
#### Input & Output file names ####
tumor_file = str_c("/media/sarthak/gender_prognosis/Tumor/",args[1],".tab")
clinical_file <- str_c("/media/sarthak/gender_prognosis/Clinical/processed/",args[1],"_pro.tab")
output_dir <- "/media/sarthak/gender_prognosis/plots/20200812_exploratory/"
#####################
tum <- read.table(file = tumor_file,  sep = '\t', header = TRUE)
tumor_samples <- colnames(tum)
rm(list="tum")
tumor_samples<- sub(".*(\\d+{4}).*$", "\\1", tumor_samples)
sample <- fread(clinical_file)
sample <- unique.data.frame(sample)
sample$ID <- sapply(strsplit(sample$submitter_id,split = "-"), "[",3)
sample$gender <- as.character(sample$gender)
sample[sample$gender=="female",]$gender <- "1"
sample[sample$gender=="male",]$gender <- "0"
sample$gender <- as.numeric(sample$gender)
sample <- sample[!is.na(sample$gender),]
sample <- sample[sample$ID %in% tumor_samples,]

gender_plot <- str_c(output_dir,args[1],"_gender.jpeg")
age_plot <- str_c(output_dir,args[1],"_age.jpeg")
jpeg(filename = gender_plot)
hist(sample$gender,main=str_c("gender_",args[1]),xlab="Gender")
dev.off()
jpeg(filename = age_plot)
hist(sample$age,main=str_c("age_",args[1]),xlab="Age(in yrs")
dev.off()
