rm(list = ls())

library(data.table)
library(stringr)
library(survival)


args <- commandArgs(TRUE)
print(args)
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/")

#### Input & Output file names ####
#tumor_file = str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Tumor/",args[1],".tab")
clinical_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Clinical/processed/",args[1],"_pro.tab")
normalized_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/data/",args[1],".tab")
output_file1 <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/",args[1],"_SPG_PG_both.tab")
output_file2 <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/analysis/prognostic/trials/",args[1],"_SPG_PG_est.tab")
#####################

##Initialising Tumor Dataset ##
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
train_rows <- sample(1:dim(data_org)[1], 0.60*dim(data_org)[1])
data_train <- data_org[train_rows, ]
data_test <- data_org[-train_rows, ]


##### Survival analysis: Cox proportional Hazard Model

#looped over all genes
data = data_train
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
rm(list = c('a','b','x','y','cox.gene','cox.gene_int'))
###################### Now preparing for LASSO ####
## for int
mean_abs1 <- quantile(abs(df_int$Interaction_coef),0.8)
length(which(abs(df_int$Interaction_coef) > mean_abs1))
output1 <- subset(df_int, subset=(abs(df_int$Interaction_coef) > mean_abs1)) 
dim(output1)
output1$FDR<- p.adjust(output1$Interaction_Pr...z.., "BH")
length(which(output1$FDR < 0.10))
output1 <- subset(output1, subset=(abs(output1$FDR) < 0.10))
## for est
mean_abs2 <- quantile(abs(df_est$gene_coef),0.8)
length(which(abs(df_est$gene_coef) > mean_abs2))
output2 <- subset(df_int, subset=(abs(df_est$gene_coef) > mean_abs2)) 
dim(output2)
output2$FDR<- p.adjust(output2$gene_Pr...z.., "BH")
length(which(output2$FDR < 0.10))
output2 <- subset(output2, subset=(abs(output2$FDR) < 0.10))

############# For both models ###########
ML <- subset(data, select=union(as.character(output1$gene),as.character(output2$gene)))
Y <- data[ , which(names(data) %in% c("PFI","PFI.time"))]
colnames(Y) <- c("status", "time")
Y$status <- as.numeric(as.character(Y$status))
library(glmnet)

# train_rows1 <- sample(1:dim(data)[1], 0.6*dim(data)[1])
# x.train <- as.matrix(ML[train_rows1, ])
# x.test <- as.matrix(ML[-train_rows1, ])
# 
# y.train <- as.matrix(Y[train_rows1, ])
# y.test <-  as.matrix(Y[-train_rows1, ])
x.train <- as.matrix(ML)
y.train <- as.matrix(Y)
fit <- glmnet(x.train, y.train, family="cox", alpha=0.1)
cv.fit <- cv.glmnet(x.train, y.train, type.measure="deviance",grouped=TRUE,alpha=0.1,family="cox",nfolds=100)
plot(fit, xvar="lambda")
plot(cv.fit, main="LASSO")
coef(cv.fit)
coef.min = coef(cv.fit, s = "lambda.min")
active.min = which(coef.min != 0)
predictors <- subset(x.train, select=active.min)
pred_cols <- colnames(predictors)
# pred <- predict(fit,x.test)
# pred
# plot(pred,y.test)
#merge.train.pred  <- merge.train[,c]
#col.num <- which(colnames(merge.train) %in% c)
#merge.train.pred <- subset(merge.train, select=col.num)
#col.num <- which(colnames(merge.train) %in% c)
# col.num <- which(colnames(data_org) %in% pred_cols)
# col.num
variable1 <- paste0('`',colnames(predictors),'`')
variable2 <- noquote(paste(variable1,collapse=" + "))
#variable2
cox.form <- as.formula(paste("Surv(PFI.time,as.numeric(PFI)) ~", variable2, sep = ""))
#cox.form
cox_prog <- coxph(cox.form,data=data_org)$coef
risk_calc <- function(dataframe,coefficient_df){
  risk = vector(length=nrow(dataframe))
  for(i in 1:length(coefficient_df)){
    risk = risk + (dataframe[names(coefficient_df)[i]]*coefficient_df[names(coefficient_df)[i]])
  }
  names(risk)= NULL
  #row.names(risk) = NULL
  return(risk)
}
cox_prog_source <- names(cox_prog)
cox_prog_source[cox_prog_source %in% output1$gene & cox_prog_source %in% output2$gene] <- 'both'
cox_prog_source[cox_prog_source %in% output1$gene] <- 'int'
cox_prog_source[cox_prog_source %in% output2$gene] <- 'est'

iters <- c('data_train','data_test')
frame_out <- data.frame('bootstrap'= args[2],'predictors' = paste0(names(cox_prog),collapse = ';'),
                        'coefficients' = paste0(cox_prog,collapse = ';'),
                        'model_source'= paste0(cox_prog_source,collapse = ';'))
for(i in 1:length(iters)){
  sub <- get(iters[i])
  sub[,"Risk_Score"] <- risk_calc(sub,cox_prog)
  dim(sub)
  sub[,"time"] <- sub$PFI.time
  sub[,"status"] <- sub$PFI
  f <- sub[sub$gender == 1,]
  m <- sub[sub$gender == 0,]
  med_sub <- median(sub$Risk_Score)
  med_f <- median(f$Risk_Score)
  med_m <- median(m$Risk_Score)
  sub$Risk_Score <- ifelse(sub$Risk_Score < med_f, 0, 1)
  f$Risk_Score <- ifelse(f$Risk_Score < med_f, 0, 1)
  m$Risk_Score <- ifelse(m$Risk_Score < med_m, 0, 1)
  fit_sub <- survdiff(Surv(PFI.time, PFI) ~ Risk_Score, data = sub)
  fit_f <- survdiff(Surv(PFI.time, PFI) ~ Risk_Score, data = f)
  fit_m <- survdiff(Surv(PFI.time, PFI) ~ Risk_Score, data = m)
  p_sub <- broom::glance(fit_sub)$p.value
  p_f <- broom::glance(fit_f)$p.value
  p_m <- broom::glance(fit_m)$p.value
  cens = nrow(sub[sub$status==1,])/nrow(sub[sub$status==0,])
  frame_in <- data.frame(cens,nrow(sub),nrow(f),nrow(m),p_sub,p_f,p_m)
  colnames(frame_in) <- c('cens_rat','all_size','f_size','m_size','all_LR','f_LR','m_LR')
  colnames(frame_in) <- paste0(gsub('data_','',iters[i]),'_',colnames(frame_in))
  frame_out <- cbind(frame_out,frame_in)
}
write.table(frame_out, file = output_file1, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)

rm(list = c('ML','Y','x.train','y.train','fit','cv.fit','coef.min','active.min','predictors','pred_cols','variable1','variable2','cox.form','cox_prog','cox_prog_source','frame_out','frame_in'))


##### for model est only

############# For both models ###########
ML <- subset(data, select=as.character(output2$gene))
Y <- data[ , which(names(data) %in% c("PFI","PFI.time"))]
colnames(Y) <- c("status", "time")
Y$status <- as.numeric(as.character(Y$status))
library(glmnet)

# train_rows1 <- sample(1:dim(data)[1], 0.6*dim(data)[1])
# x.train <- as.matrix(ML[train_rows1, ])
# x.test <- as.matrix(ML[-train_rows1, ])
# 
# y.train <- as.matrix(Y[train_rows1, ])
# y.test <-  as.matrix(Y[-train_rows1, ])
x.train <- as.matrix(ML)
y.train <- as.matrix(Y)
fit <- glmnet(x.train, y.train, family="cox", alpha=0.1)
cv.fit <- cv.glmnet(x.train, y.train, type.measure="deviance",grouped=TRUE,alpha=0.1,family="cox",nfolds=100)
#plot(fit, xvar="lambda")
#plot(cv.fit, main="LASSO")
coef(cv.fit)
coef.min = coef(cv.fit, s = "lambda.min")
active.min = which(coef.min != 0)
predictors <- subset(x.train, select=active.min)
pred_cols <- colnames(predictors)
# pred <- predict(fit,x.test)
# pred
# plot(pred,y.test)
#merge.train.pred  <- merge.train[,c]
#col.num <- which(colnames(merge.train) %in% c)
#merge.train.pred <- subset(merge.train, select=col.num)
#col.num <- which(colnames(merge.train) %in% c)
# col.num <- which(colnames(data_org) %in% pred_cols)
# col.num
variable1 <- paste0('`',colnames(predictors),'`')
variable2 <- noquote(paste(variable1,collapse=" + "))
#variable2
cox.form <- as.formula(paste("Surv(PFI.time,as.numeric(PFI)) ~", variable2, sep = ""))
#cox.form
cox_prog <- coxph(cox.form,data=data_org)$coef
risk_calc <- function(dataframe,coefficient_df){
  risk = vector(length=nrow(dataframe))
  for(i in 1:length(coefficient_df)){
    risk = risk + (dataframe[names(coefficient_df)[i]]*coefficient_df[names(coefficient_df)[i]])
  }
  names(risk)= NULL
  #row.names(risk) = NULL
  return(risk)
}
cox_prog_source <- 'est_model_only'

iters <- c('data_train','data_test')
frame_out <- data.frame('bootstrap'= args[2],'predictors' = paste0(names(cox_prog),collapse = ';'),
                        'coefficients' = paste0(cox_prog,collapse = ';'),
                        'model_source'= paste0(cox_prog_source,collapse = ';'))
for(i in 1:length(iters)){
  sub <- get(iters[i])
  sub[,"Risk_Score"] <- risk_calc(sub,cox_prog)
  dim(sub)
  sub[,"time"] <- sub$PFI.time
  sub[,"status"] <- sub$PFI
  f <- sub[sub$gender == 1,]
  m <- sub[sub$gender == 0,]
  med_sub <- median(sub$Risk_Score)
  med_f <- median(f$Risk_Score)
  med_m <- median(m$Risk_Score)
  sub$Risk_Score <- ifelse(sub$Risk_Score < med_f, 0, 1)
  f$Risk_Score <- ifelse(f$Risk_Score < med_f, 0, 1)
  m$Risk_Score <- ifelse(m$Risk_Score < med_m, 0, 1)
  fit_sub <- survdiff(Surv(PFI.time, PFI) ~ Risk_Score, data = sub)
  fit_f <- survdiff(Surv(PFI.time, PFI) ~ Risk_Score, data = f)
  fit_m <- survdiff(Surv(PFI.time, PFI) ~ Risk_Score, data = m)
  p_sub <- broom::glance(fit_sub)$p.value
  p_f <- broom::glance(fit_f)$p.value
  p_m <- broom::glance(fit_m)$p.value
  cens = nrow(sub[sub$status==1,])/nrow(sub[sub$status==0,])
  frame_in <- data.frame(cens,nrow(sub),nrow(f),nrow(m),p_sub,p_f,p_m)
  colnames(frame_in) <- c('cens_rat','all_size','f_size','m_size','all_LR','f_LR','m_LR')
  colnames(frame_in) <- paste0(gsub('data_','',iters[i]),'_',colnames(frame_in))
  frame_out <- cbind(frame_out,frame_in)
}
write.table(frame_out, file = output_file2, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
