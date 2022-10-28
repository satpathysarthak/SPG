### Validating ressssssults on CGGA data 
rm(list = ls())

library(data.table)
library(stringr)
# #### Input & Output file names ####
# tumor_file = '/home/workstation/Documents/Sarthak/gender_prognosis/validation/cgga_tumor.txt'
# clinical_file = '/home/workstation/Documents/Sarthak/gender_prognosis/validation/cgga_clinical_GBM.txt'
# output_file = '/home/workstation/Documents/Sarthak/gender_prognosis/validation/cgga_data_gbm.txt'
# 
# #####initialising Sample Dataset#####
# sample <- fread(clinical_file)
# sample <- unique.data.frame(sample)
# #sample$ID <- sapply(strsplit(sample$submitter_id,split = "-"), "[",3)
# sample$Gender <- as.character(sample$Gender)
# sample[sample$Gender=="Female",]$Gender <- "1"
# sample[sample$Gender=="Male",]$Gender <- "0"
# sample$Gender <- as.numeric(sample$Gender)
# sample <- sample[!is.na(sample$Gender),]
# 
# ##Initialising Tumor Dataset ##
# tum <- read.table(file = tumor_file,  sep = '\t', header = TRUE)
# rownames(tum) <- tum$Gene_Name
# tum$Gene_Name <- NULL
# t_row <- rownames(tum)
# t_col <- colnames(tum)
# tum <- transpose(tum)
# rownames(tum) <- t_col
# colnames(tum) <- t_row
# #tum$CGGA_ID <- rownames(tum)
# tum <- tum[rownames(tum) %in% sample$CGGA_ID,]
# normalize <- function(df, cols) {
#   result <- df # make a copy of the input data frame
#   
#   for (j in cols) { # each specified col
#     m <- mean(df[,j]) # column mean
#     std <- sd(df[,j]) # column (sample) sd
#     
#     for (i in 1:nrow(result)) { # each row of cur col
#       result[i,j] <- (result[i,j] - m) / std
#     }
#   }
#   return(result)
# }
# tum <- normalize(tum,colnames(tum))
# tum$CGGA_ID <- rownames(tum)
# ##### Data Merging #####
# data <- merge(sample, tum, by="CGGA_ID", all = F)
# data <- data.frame(data)
# write.table(data, output_file, append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')

## Calculating Risk Scores

# required files
input_file <- '/home/workstation/Documents/Sarthak/gender_prognosis/validation/cgga_data_gbm.txt'
model_boot_both_file <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/GBM_SPG_PG_both.tab' 
model_boot_est_file <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/GBM_SPG_PG_est.tab' 
map_genes <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/GBM_model_genes_mapped_97.txt'

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
model_no <- 97
cox_prog_both <- as.numeric(strsplit(model_both$coefficients[model_both$bootstrap==model_no],';')[[1]])
names(cox_prog_both) <- strsplit(model_both$predictors[model_both$bootstrap==model_no],';')[[1]]
cox_prog_est <- as.numeric(strsplit(model_est$coefficients[model_est$bootstrap==model_no],';')[[1]])
names(cox_prog_est) <- strsplit(model_est$predictors[model_est$bootstrap==model_no],';')[[1]]

df_both_genes <- data.frame('gene' = names(cox_prog_both),'coef'=cox_prog_both)
df_est_genes <- data.frame('gene' = names(cox_prog_est), 'coef'=cox_prog_est)

df_both_genes <- merge(df_both_genes, map, by='gene', all.y=F)
df_est_genes <- merge(df_est_genes, map, by='gene', all.y=F)

df_both_genes <- df_both_genes[df_both_genes$gene_name != '',]
df_est_genes <- df_est_genes[df_est_genes$gene_name != '',]

cox_prog_both <- df_both_genes$coef
names(cox_prog_both) <- df_both_genes$gene_name
cox_prog_est <- df_est_genes$coef
names(cox_prog_est) <- df_est_genes$gene_name

data = data_org
data$Risk_Score_both <- risk_calc(data,cox_prog_both)
data$Risk_Score_est <- risk_calc(data,cox_prog_est)
foo <- data.frame()
for(i in 1:8){
  auc1 <- survivalROC(Stime=data$OS, status=data$Censor..alive.0..dead.1.,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")$AUC
  auc2 <- survivalROC(Stime=data$OS, status=data$Censor..alive.0..dead.1.,marker=data$Risk_Score_est,predict.time=(i*365),method="KM")$AUC
  #print(str_c(i,';',auc1,';',auc2))
  foo <- rbind(foo,data.frame('year'=i,'auc1'=auc1,'auc2'=auc2))
}
i=3
auc1 <- survivalROC(Stime=data$OS, status=data$Censor..alive.0..dead.1.,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
auc2 <- survivalROC(Stime=data$OS, status=data$Censor..alive.0..dead.1.,marker=data$Risk_Score_est,predict.time=(i*365),method="KM")
plot(auc1$FP,auc1$TP,type="l",xlim=c(0,1),ylim=c(0,1),xlab="FP",ylab="TP",col=c("#d92128"),lwd=2)
lines(auc2$FP,auc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#0c8c00"),lwd=2)
legend('bottomright', legend=c("Risk Score Est+Int", "Risk Score Est"),col=c("#d92128","#0c8c00"), cex=0.5,lty=1,lwd=2)
abline(0,1,lwd=2)

data = data_org
data$Risk_Score_both <- risk_calc(data,cox_prog_both)
data_m <- subset(data,Gender==0)
data_f <- subset(data,Gender==1)
#data$Risk_Score_est <- risk_calc(data,cox_prog_est)

# for(i in 1:16){
#   auc1 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")$AUC
#   auc2 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_est,predict.time=(i*365),method="KM")$AUC
#   print(str_c(i,';',auc1,';',auc2))
# }
i=3
data=data_m
auc1 <- survivalROC(Stime=data$OS, status=data$Censor..alive.0..dead.1.,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
data=data_f
auc2 <- survivalROC(Stime=data$OS, status=data$Censor..alive.0..dead.1.,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
plot(auc1$FP,auc1$TP,type="l",xlim=c(0,1),ylim=c(0,1),xlab="FP",ylab="TP",col=c("#215ed9"),lwd=2)
lines(auc2$FP,auc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#d9219c"),lwd=2)
legend('bottomright', legend=c("Male", "Female"),col=c("#215ed9","#d9219c"), cex=0.5,lty=1,lwd=2)
abline(0,1,lwd=2)

df <- data.frame()
for(i in 1:8){
  data=data_m
  auc1 <- survivalROC(Stime=data$OS, status=data$Censor..alive.0..dead.1.,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
  data=data_f
  auc2 <- survivalROC(Stime=data$OS, status=data$Censor..alive.0..dead.1.,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
  df <- rbind(df,data.frame('year'=i,'auc_m'=auc1$AUC,'auc_f'=auc2$AUC))
}
