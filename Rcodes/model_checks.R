# Evaluates the models generated from bootstrap method
rm(list = ls())

library(data.table)
library(stringr)
library(dplyr)
library(survival)
library(survivalROC)

args <- commandArgs(TRUE)
print(args)
setwd('/home/sarthak/Documents/projects/MKlab/')

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
#tumor_file = str_c("/home/sarthak/Documents/projects/MKlab/Tumor/",args[1],".tab")
clinical_file <- str_c("/home/sarthak/Documents/projects/MKlab/Clinical/processed/",args[1],"_pro.tab")
normalized_file <- str_c("/home/sarthak/Documents/projects/MKlab/analysis/prognostic/data/",args[1],".tab")
model_boot_both_file <- '/home/sarthak/Documents/projects/MKlab/analysis/prognostic/trials/HNSC_SPG_PG_both.tab' 
model_boot_est_file <- '/home/sarthak/Documents/projects/MKlab/analysis/prognostic/trials/HNSC_SPG_PG_est.tab' 
#output_file <- str_c("/home/sarthak/Documents/projects/MKlab/analysis/prognostic/trials/",args[1],"_SPG.tab")
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
model_no

model_no <- 33
cox_prog_both <- as.numeric(strsplit(model_both$coefficients[model_both$bootstrap==model_no],';')[[1]])
names(cox_prog_both) <- strsplit(model_both$predictors[model_both$bootstrap==model_no],';')[[1]]
cox_prog_est <- as.numeric(strsplit(model_est$coefficients[model_est$bootstrap==model_no],';')[[1]])
names(cox_prog_est) <- strsplit(model_est$predictors[model_est$bootstrap==model_no],';')[[1]]


data = data_org
data$Risk_Score_both <- risk_calc(data,cox_prog_both)
data$Risk_Score_est <- risk_calc(data,cox_prog_est)
data$Risk_Score_both <- unlist(data$Risk_Score_both)
data$Risk_Score_est <- unlist(data$Risk_Score_est)


for(i in 1:8){
  auc1 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")$AUC
  auc2 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_est,predict.time=(i*365),method="KM")$AUC
  print(str_c(i,';',auc1,';',auc2))
}
i=4
auc1 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
auc2 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_est,predict.time=(i*365),method="KM")
plot(auc1$FP,auc1$TP,type="l",xlim=c(0,1),ylim=c(0,1),xlab="FP",ylab="TP",col=c("#d92128"),lwd=2,cex.lab=1.25, cex.axis=1.25)
lines(auc2$FP,auc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#0c8c00"),lwd=2)
legend('bottomright', legend=c("Risk Score Est+Int", "Risk Score Est"),col=c("#d92128","#0c8c00"), cex=1.5,lty=1,lwd=2)
abline(0,1,lwd=2)


model_no <- 1
foo <- merge(model_both,model_est,by='bootstrap')
foo2 <- data.frame()
foo3 <- foo[foo$train_all_LR.x< 0.05 & foo$test_all_LR.x,]
for(model_no in foo3$bootstrap){
  cox_prog_both <- as.numeric(strsplit(model_both$coefficients[model_both$bootstrap==model_no],';')[[1]])
  names(cox_prog_both) <- strsplit(model_both$predictors[model_both$bootstrap==model_no],';')[[1]]
  cox_prog_est <- as.numeric(strsplit(model_est$coefficients[model_est$bootstrap==model_no],';')[[1]])
  names(cox_prog_est) <- strsplit(model_est$predictors[model_est$bootstrap==model_no],';')[[1]]
  
  
  data = data_org
  data$Risk_Score_both <- risk_calc(data,cox_prog_both)
  data$Risk_Score_est <- risk_calc(data,cox_prog_est)
  
  # for(i in 1:16){
  #   auc1 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")$AUC
  #   auc2 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_est,predict.time=(i*365),method="KM")$AUC
  #   print(str_c(i,';',auc1,';',auc2))
  # }
  i=6
  auc1 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
  auc2 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_est,predict.time=(i*365),method="KM")
  #print(str_c(model_no,';',auc1$AUC,';',auc2$AUC))
  foo2 <- rbind(foo2,data.frame(model_no,auc1$AUC,auc2$AUC))
}
foo <- merge(model_both,model_est,by='bootstrap')
foo <- foo[foo$train_all_LR.x<0.05&foo$test_all_LR.y<0.05,]
foo$predictors_no.x <- NA
for(i in 1:nrow(foo)){
  foo$predictors_no.x[i]<-length(sapply(strsplit(foo$predictors.x[i], ';'),'['))
}
foo$predictors_no.y <- NA
for(i in 1:nrow(foo)){
  foo$predictors_no.y[i]<-length(sapply(strsplit(foo$predictors.y[i], ';'),'['))
}
library(ggplot2)
ggplot(aes(predictors_no.x), data = foo) + 
  geom_histogram(color="black", fill="white", bins = 10)
gene_list <- paste0(foo$predictors.x, collapse = ';') 
gene_list <- sapply(strsplit(gene_list, ';'),'[')
colnames(gene_list) <- 'gene'
gene_list <- as.data.frame(gene_list)
gene_list <- unique.data.frame(gene_list)
gene_list$gene <- as.character(gene_list$gene)
gene_count <- function(gene, frame){
  return(sum(str_count(frame,gene)))
}
gene_count(gene_list$gene[1],foo$predictors.x)
gene_list$gene_both <- NA
gene_list$gene_est <- NA

for(i in 1:nrow(gene_list)){
  gene_list$gene_both[i] <- gene_count(gene_list$gene[i],foo$predictors.x)
  gene_list$gene_est[i] <- gene_count(gene_list$gene[i],foo$predictors.y)
}
gene_list$total <- gene_list$gene_both + gene_list$gene_est
gene_map <- fread('/home/sarthak/Documents/projects/MKlab/analysis/prognostic/trials/gene-list_mapped')
gene_map <- gene_map[,c(1,8,9)]
long <- gene_list[,c(1:4)]
long <- merge(long,gene_map,by='gene')
for(i in 1:nrow(long)){
  if(long$gene_name[i]==''){long$gene_name[i]=long$gene_ensembl[i]}
}
long$gene <- long$gene_name
long$gene_name <- NULL
long$gene_ensembl <- NULL
long <- reshape2::melt(data.frame(long),id=c('gene','total'))
long %>%
  arrange(desc(gene))%>%
  arrange(desc(total))%>%
  head(40)%>%
  ggplot(aes(fill=variable, y=value, x=reorder(gene,-total))) + 
  geom_bar(position="stack", stat="identity")
## counting by genes: est, int both
frame1 <- foo$predictors.x
frame2 <- foo$model_source.x
df_gene <- data.frame()
for(j in 1:length(gene_list$gene)){
  gene <- gene_list$gene[j]
  arr_gene <- vector()
  for(i in 1:nrow(foo)){
    if(grepl(gene,frame1[i]) == T){
      arr1 <- strsplit(frame1[i], ';')[[1]]
      arr2 <- strsplit(frame2[i], ';')[[1]]
      source <- arr2[grep(gene,arr1)]
      arr_gene <- append(arr_gene,source)
    }
  }
  df_gene <- rbind(df_gene,data.frame('gene'=gene, 'est'=length(grep('est',arr_gene)),
             'int'=length(grep('int',arr_gene)),
             'both'=length(grep('both',arr_gene))))
}
gene_list <- merge(gene_list,df_gene, by='gene')

gene_map <- fread('/home/sarthak/Documents/projects/MKlab/analysis/prognostic/trials/gene-list_mapped')
gene_map <- gene_map[,c(1,8,9)]
long <- gene_list
long$total= NULL
long$gene_est = NULL

long <- merge(long,gene_map,by='gene')
for(i in 1:nrow(long)){
  if(long$gene_name[i]==''){long$gene_name[i]=long$gene_ensembl[i]}
}
long$gene <- long$gene_name
long$gene_name <- NULL
long$gene_ensembl <- NULL


long <- reshape2::melt(data.frame(long),id=c('gene','gene_both'))
long %>%
  arrange(desc(gene))%>%
  arrange(desc(gene_both))%>%
  head(60)%>%
  ggplot(aes(fill=variable, y=value, x=reorder(gene,-gene_both))) + 
  geom_bar(position="stack", stat="identity")

top_gene <- long %>%
  arrange(desc(gene))%>%
  arrange(desc(total))%>%
  head(40) 
top_gene <- unique(top_gene$gene)
frame1 <- foo$predictors.x
foo$overlap <- NA
for(i in 1:nrow(foo)){
  arr1 <- strsplit(frame1[i], ';')[[1]]
  foo$overlap[i] <- length(intersect(arr1,top_gene))
  #arr_gene <- append(arr_gene,source)
}
foo$overlap <- foo$overlap/foo$predictors_no.x
  
#####
model_no <- 97
cox_prog_both <- as.numeric(strsplit(model_both$coefficients[model_both$bootstrap==model_no],';')[[1]])
names(cox_prog_both) <- strsplit(model_both$predictors[model_both$bootstrap==model_no],';')[[1]]
#cox_prog_est <- as.numeric(strsplit(model_est$coefficients[model_est$bootstrap==model_no],';')[[1]])
#names(cox_prog_est) <- strsplit(model_est$predictors[model_est$bootstrap==model_no],';')[[1]]


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
i=6
data=data_m
auc3 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
data=data_f
auc4 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
lines(auc3$FP,auc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#92C5DE"),lwd=2)
lines(auc4$FP,auc4$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#F1B6DA"),lwd=2)
#legend('bottomright', legend=c("Male", "Female"),col=c("#92C5DE","#F1B6DA"), cex=0.5,lty=1,lwd=2)
legend('bottomright', legend=c("Risk Score Est+Int", "Risk Score Est","Risk Male Est+Int", "Risk Female Est+Int"),col=c("#d92128","#0c8c00","#92C5DE","#F1B6DA"),cex=0.75,lty=1,lwd=2)
text()
abline(0,1,lwd=2)

df <- data.frame()
for(i in 1:8){
  data=data_m
  auc1 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
  data=data_f
  auc2 <- survivalROC(Stime=data$PFI.time, status=data$PFI,marker=data$Risk_Score_both,predict.time=(i*365),method="KM")
  df <- rbind(df,data.frame('year'=i,'auc_m'=auc1$AUC,'auc_f'=auc2$AUC))
}
ggplot()
