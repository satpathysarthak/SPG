## Script to perform Survival Analysis with Tumor data using 2 models: Male Female
rm(list = ls())
args<-commandArgs(TRUE)
print(args)
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/")
library(data.table)
library(stringr)
#### Input & Output file names ####
tumor_file = str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Tumor/",args[1],".tab")
clinical_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Clinical/processed/",args[1],"_pro.tab")
fdr_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore_name/",args[1],"_FDR_int.name")
output_file_male <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/bivariate_gender/",args[1],"_male.tab")
output_file_female <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/bivariate_gender/",args[1],"_female.tab")
#####################

##Initialising Tumor Dataset ##
tum <- read.table(file = tumor_file,  sep = '\t', header = TRUE)
fdr <- fread(fdr_file,header = T)
tum <- tum[fdr$gene,]
# length_dim <- nrow(tum)
# width_dim <- ncol(tum)
# median_tum <- vector()
# for (i in 1:length_dim)#Till the last gene entry
# {
#   median_tum[i] <- median(as.matrix(tum[i,]))#loop is run only till the last patient entry
#   #if the covariates are not mentioned the loop is not functional after the 1st iteration
#   if(i%%5000==0){
#     print(i)
#   }
# }
# tum$median <- median_tum
# tum <- subset(tum, tum$median >1)#filtering the ones with more than 1FPKM
# tum$median = NULL
# median_tum <- NULL
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
#a =tum[,c(1:3)]
#b = normalize(a,colnames(a))
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
##### Survival analysis: Cox proportional Hazard Model
library(survival)
#looped over all genes
rm(list = "tum")
gen_i <- grep("^gender", colnames(data))
age_i <- grep("^age", colnames(data))
low <- ncol(sample) + 1
high <- ncol(data)
df_male <- data.frame()
df_female <- data.frame()
data_male <- data[data$gender == 0,]
data_female <- data[data$gender == 1,]
for (i in low:high)#: range of genes in the data(merged dataframe)
{
  #Block for taking gene(individually), age and gender
  #4th column is age and 5th is gender
  cox.gene <- coxph(Surv(PFI.time,PFI)~data_male[,i]+data_male[,"age"], data=data_male)
  a <- data.frame(coef(summary(cox.gene))[,c(1:5)])#sores only the coeff, exp, z , pvalue from the model
  b <- cbind(colnames(data_male)[i], a[1,], colnames(data_male)[age_i], a[2,])# melts the 3 rows to one entry by cbind
  rownames(b) <- colnames(data_male)[i]
  colnames(b)[1] <- 'gene'
  colnames(b)[2:6] <- str_c('gene_',colnames(b)[2:6])
  colnames(b)[7] <- 'age'
  colnames(b)[8:12] <- str_c('age_',colnames(b)[8:12])
  #storing the output in a .tab file
  #Set the working directory or provide the path
  df_male <- rbind(df_male,b)
  
  #Similarly, the block that incorporates the interaction term age*gender
  cox.gene_int <- coxph(Surv(PFI.time,PFI)~data_female[,i]+data_female[,"age"], data = data_female)
  x <- data.frame(coef(summary(cox.gene_int))[,c(1:5)])
  y <- cbind(colnames(data_female)[i], x[1,], colnames(data_female)[age_i], x[2,])
  rownames(y) <- colnames(data)[i]
  colnames(y)[1] <- 'gene'
  colnames(y)[2:6] <- str_c('gene_',colnames(y)[2:6])
  colnames(y)[7] <- 'age'
  colnames(y)[8:12] <- str_c('age_',colnames(y)[8:12])
  df_female <- rbind(df_female,y)
  
}
write.table(df_male, output_file_male, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
write.table(df_female, output_file_female, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
warnings()
rm(list = ls())