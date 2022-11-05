## Script to generate plots for data exploration: ditribution for cencored and gender data`
rm(list = ls())
args<-commandArgs(TRUE)
print(args)
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/")
library(data.table)
library(stringr)
#### Input & Output file names ####
tumor_file = str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Tumor/",args[1],".tab")
clinical_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Clinical/processed/",args[1],"_pro.tab")
output_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/table"
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
num_PFI_0 <- nrow(sample[sample$PFI==0,])
num_PFI_1 <- nrow(sample[sample$PFI==1,])
num_gen_0 <- nrow(sample[sample$gender==0,])
num_gen_1 <- nrow(sample[sample$gender==1,])
out <- data.frame("cancer"=args[1],"patients"=length(unique(sample$submitter_id)),
                  "PFI_0"=num_PFI_0,"PFI_1"=num_PFI_1,
                  "male"=num_gen_0,"male"=num_gen_1)
write.table(out, file = str_c(output_dir,"/20201027_list.txt"),
            append = T, col.names = F, sep = "\t", row.names = F) 
