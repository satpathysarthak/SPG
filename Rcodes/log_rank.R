#### Calculating Log rank for all FDR genes

args <- commandArgs(TRUE)
library("survival")
library("survminer")
library(data.table)
library(stringr)
#### Input & Output file names ####
tumor_file = str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Tumor/",args[1],".tab")
clinical_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Clinical/processed/",args[1],"_pro.tab")
fdr_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore_name/",args[1],"_FDR_int.name")
out_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/analysis/log_rank/", args[1],"_log_rank_int.name")
#####################

##Initialising Tumor Dataset ##

tum <- read.table(file = tumor_file,  sep = '\t', header = TRUE)
fdr <- fread(fdr_file,header = T)
tum <- tum[fdr$gene,]
t_row <- rownames(tum)
t_col <- colnames(tum)
tum <- transpose(tum)
rownames(tum) <- t_col
colnames(tum) <- t_row
#tum <- as.data.frame(tum)
normalize <- function(df, cols) {
  result <- df # make a copy of the input data frame
  
  for (j in cols) { # each specified col
    m <- mean(df[,j]) # column mean
    std <- sd(df[,j]) # column (sample) sd
    
    for (i in 1:nrow(result)) { # each row of cur col
      result[i,j] <- (result[i,j] - m) / std
    }
  }
  return(result)
}
tum <- normalize(tum,colnames(tum))
#tum <- normalize(tum[,c(1:3)],colnames(tum)[c(1:3)])
tum$ID <- rownames(tum)
tum$ID <- sub(".*(\\d+{4}).*$", "\\1", tum$ID)

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
data <- merge(sample, tum, by="ID", all = F)
data <- data.frame(data)

male_sub <- data[data$gender==0,]
female_sub <- data[data$gender==1,]

gen_i <- grep("^gender", colnames(data))
age_i <- grep("^age", colnames(data))
low <- ncol(sample) + 1
high <- ncol(data)
data_log_rank <- data.frame()
for (j in low:high){
  med_male <- median(male_sub[,j])
  med_female <- median(female_sub[,j])
  male_sub[,j]<-ifelse(male_sub[,j]< med_male, 0, 1)
  female_sub[,j]<-ifelse(female_sub[,j]< med_female, 0, 1)
  frame_gene <- data.frame("gene"=colnames(data)[j])  #gene name
  male_fit <- survdiff(Surv(PFI.time, PFI) ~ male_sub[,j], data = male_sub)
  male_p <- broom::glance(male_fit)$p.value
  frame_male <- data.frame(male_fit$n[1],male_fit$n[2],male_fit$obs[1],male_fit$obs[2],male_fit$exp[1],male_fit$exp[2],male_p)
  rownames(frame_male) = NULL
  frame_male$male_fit.diffs.1. = (male_fit$obs[1]-male_fit$exp[1])**2/male_fit$exp[1]
  frame_male$male_fit.diffs.2. = (male_fit$obs[2]-male_fit$exp[2])**2/male_fit$exp[2]
  female_fit <- survdiff(Surv(PFI.time, PFI) ~ female_sub[,j], data = female_sub)
  female_p <- broom::glance(female_fit)$p.value
  frame_female <- data.frame(female_fit$n[1],female_fit$n[2],female_fit$obs[1],female_fit$obs[2],female_fit$exp[1],female_fit$exp[2],female_p)
  rownames(frame_female) = NULL
  frame_female$female_fit.diffs.1. = (female_fit$obs[1]-female_fit$exp[1])**2/female_fit$exp[1]
  frame_female$female_fit.diffs.2. = (female_fit$obs[2]-female_fit$exp[2])**2/female_fit$exp[2]
  data_lr_gene <- cbind(frame_gene,frame_male,frame_female)
  data_log_rank <- rbind(data_log_rank,data_lr_gene)
  rm(list=c("data_lr_gene","frame_male","frame_female","frame_gene","male_fit","female_fit","male_p","female_p"))
}
write.table(data_log_rank, out_file, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
warnings()