## Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated
rm(list = ls())
args<-commandArgs(TRUE)
print(args)
setwd("/media/sarthak/gender_prognosis/")
library(data.table)
library(stringr)
#### Input & Output file names ####
tumor_file = str_c("/media/sarthak/gender_prognosis/Tumor/",args[1],".tab")
clinical_file <- str_c("/media/sarthak/gender_prognosis/Clinical/processed/",args[1],"_pro.tab")
output_file_est <- str_c("/media/sarthak/gender_prognosis/cox_regression/gender/",args[1],"_est.tab")
output_file_int <- str_c("/media/sarthak/gender_prognosis/cox_regression/gender/",args[1],"_int.tab")
#####################

##Initialising Tumor Dataset ##
tum <- read.table(file = tumor_file,  sep = '\t', header = TRUE)
length_dim <- nrow(tum)
width_dim <- ncol(tum)
median_tum <- vector()
for (i in 1:length_dim)#Till the last gene entry
{
   median_tum[i] <- median(as.matrix(tum[i,]))#loop is run only till the last patient entry
   #if the covariates are not mentioned the loop is not functional after the 1st iteration
   if(i%%5000==0){
     print(i)
   }
}
tum$median <- median_tum
tum <- subset(tum, tum$median >1)#filtering the ones with more than 1FPKM
tum$median = NULL
median_tum <- NULL
t_row <- rownames(tum)
t_col <- colnames(tum)
tum <- transpose(tum)
rownames(tum) <- t_col
colnames(tum) <- t_row
#tum <- as.data.frame(tum)
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
##### Survival analysis: Cox proportional Hazard Model
library(survival)
#looped over all genes
rm(list = "tum")
gen_i <- grep("^gender", colnames(data))
age_i <- grep("^age", colnames(data))
low <- ncol(sample) + 1
high <- ncol(data)
df_est <- data.frame()
df_int <- data.frame()
for (i in low:high)#: range of genes in the data(merged dataframe)
{
  #Block for taking gene(individually), age and gender
  #4th column is age and 5th is gender
  cox.gene <- coxph(Surv(PFI.time,PFI)~data[,i]+data[,"age"]+data[,"gender"], data=data)
  a <- data.frame(coef(summary(cox.gene))[,c(1:5)])#sores only the coeff, exp, z , pvalue from the model
  b <- cbind(colnames(data)[i], a[1,], colnames(data)[age_i], a[2,], colnames(data)[gen_i], a[3,])# melts the 3 rows to one entry by cbind
  rownames(b) <- colnames(data)[i]
  colnames(b)[1] <- 'gene'
  colnames(b)[2:6] <- str_c('gene_',colnames(b)[2:6])
  colnames(b)[7] <- 'age'
  colnames(b)[8:12] <- str_c('age_',colnames(b)[8:12])
  colnames(b)[13] <- 'gender'
  colnames(b)[14:18] <- str_c('gender_',colnames(b)[14:18])
  #storing the output in a .tab file
  #Set the working directory or provide the path
  df_est <- rbind(df_est,b)
  
  #Similarly, the block that incorporates the interaction term age*gender
  cox.gene_int <- coxph(Surv(PFI.time,PFI)~data[,i]+data[,"age"]+data[,"gender"]+data[,i]*data[,"gender"], data = data)
  x <- data.frame(coef(summary(cox.gene_int))[,c(1:5)])
  y <- cbind(colnames(data)[i], x[1,], colnames(data)[age_i], x[2,], colnames(data)[gen_i], x[3,], "Interaction", x[4,])
  #rownames(y) <- colnames(data)[i]
  colnames(y)[1] <- 'gene'
  colnames(y)[2:6] <- str_c('gene_',colnames(y)[2:6])
  colnames(y)[7] <- 'age'
  colnames(y)[8:12] <- str_c('age_',colnames(y)[8:12])
  colnames(y)[13] <- 'gender'
  colnames(y)[14:18] <- str_c('gender_',colnames(y)[14:18])
  colnames(y)[19] <- 'Interaction'
  colnames(y)[20:24] <- str_c('Interaction_',colnames(y)[20:24])
  df_int <- rbind(df_int,y)
  
}
write.table(df_est, output_file_est, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
write.table(df_int, output_file_int, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
warnings()
rm(list = ls())