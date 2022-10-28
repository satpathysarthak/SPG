### AUC for HNSC new data
library(data.table)
library(survivalROC)
setwd('/home/sarthak/Documents/projects/MKlab')
data <- fread('output/data_HNSC.tab', header = T)


# required files
input_file <- '/home/sarthak/Documents/projects/MKlab/validation/HNSC/output/data_HNSC.tab'
model_boot_both_file <- '/home/sarthak/Documents/projects/MKlab/analysis/prognostic/trials/HNSC_1_SPG_PG_both.tab' 
model_boot_est_file <- '/home/sarthak/Documents/projects/MKlab/analysis/prognostic/trials/HNSC_1_SPG_PG_est.tab' 
map_genes <- '/home/sarthak/Documents/projects/MKlab/analysis/prognostic/trials/HNSC_model_genes_mapped.txt'

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
model_no <- 33
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

data = data_org
data$Risk_Score_both <- risk_calc(data,cox_prog_both)
data$Risk_Score_est <- risk_calc(data,cox_prog_est)
foo <- data.frame()
for(i in 1:8){
  auc1 <- survivalROC(Stime=data$pfs, status=as.numeric(data$pfs_event),marker=data$Risk_Score_both,predict.time=(i*365),method="KM")$AUC
  auc2 <- survivalROC(Stime=data$pfs, status=as.numeric(data$pfs_event),marker=data$Risk_Score_est,predict.time=(i*365),method="KM")$AUC
  #print(str_c(i,';',auc1,';',auc2))
  foo <- rbind(foo,data.frame('year'=i,'auc1'=auc1,'auc2'=auc2))
}
i=4
auc1 <- survivalROC(Stime=data$pfs, status=as.numeric(data$pfs_event),marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
auc2 <- survivalROC(Stime=data$pfs, status=as.numeric(data$pfs_event),marker=data$Risk_Score_est,predict.time=(i*365),method="KM")
plot(auc1$FP,auc1$TP,type="l",xlim=c(0,1),ylim=c(0,1),xlab="FP",ylab="TP",col=c("#d92128"),lwd=2)
lines(auc2$FP,auc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#0c8c00"),lwd=2)
legend('bottomright', legend=c("Risk Score Est+Int", "Risk Score Est"),col=c("#d92128","#0c8c00"), cex=0.5,lty=1,lwd=2)
abline(0,1,lwd=2)

data = data_org
data$Risk_Score_both <- risk_calc(data,cox_prog_both)
data_m <- subset(data,gender==0)
data_f <- subset(data,gender==1)
#data$Risk_Score_est <- risk_calc(data,cox_prog_est)

# for(i in 1:16){
#   auc1 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")$AUC
#   auc2 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_est,predict.time=(i*365),method="KM")$AUC
#   print(str_c(i,';',auc1,';',auc2))
# }
i=4
data=data_m
auc1 <- survivalROC(Stime=data$pfs, status=as.numeric(data$pfs_event),marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
data=data_f
auc2 <- survivalROC(Stime=data$pfs, status=as.numeric(data$pfs_event),marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
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

