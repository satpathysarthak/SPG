setwd("/home/workstation/Documents/Sarthak/gender_prognosis")
args <- commandArgs(TRUE)
list_genes = args[2]
#args<-"BLCA"
#print(args[1])
print(list_genes)
load("/home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/idmap_fun.RData")
library(plyr)
library(dplyr)
library(readr)
library(stringr)
library('survival')
library('survminer')
library('grid')
library('gridExtra')
library(data.table)

#### Input & Output file names ####
tumor_file = str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Tumor/",args[1],".tab")
clinical_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/Clinical/processed/",args[1],"_pro.tab")
dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210427_KM/"
fdr_file <- str_c("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore_name/",args[1],"_FDR_int.name")

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

############### Plot variables change genes #########
j = grep(list_genes[1],colnames(data))

med_male <- median(male_sub[,j])
med_female <- median(female_sub[,j])

print(med_male)
print(med_female)
male_sub[,j]<-ifelse(male_sub[,j]< med_male, 0, 1)
female_sub[,j]<-ifelse(female_sub[,j]< med_female, 0, 1)
# PFI_time_col <- grep("^PFI.time$", colnames(data))
# PFI_col <- grep("^PFI$", colnames(data))
fit_m <- survfit( Surv(time = PFI.time,
                       event = PFI)~ male_sub[,j],
                  data = male_sub)
fit_f <- survfit( Surv(time = PFI.time,
                       event = PFI)~ female_sub[,j],
                  data = female_sub)
tit_male <- str_c(args[1],"_Male_", idmap(colnames(male_sub)[j]))
tit_female <- str_c(args[1],"_Female_", idmap(colnames(female_sub)[j]))
tit_plot <- str_c(dir, args[1], colnames(data)[j],".png")
p1m <- ggsurvplot(fit_m, data=male_sub,palette = c("#404040", "#ca0020"),
                  ggtheme = theme_bw(),pval = TRUE,
                  pval.method=TRUE,conf.int = FALSE,
                  font.main=c(20,"bold"),font.x=c(16,"bold"),font.y=c(16,"bold"),
                  font.tickslab=c(14,"plain"), title = tit_male, legend.labs = c("Low", "High"))
p1f <- ggsurvplot(fit_f, data=female_sub,palette = c("#404040", "#ca0020"),
                  ggtheme = theme_bw(),pval = TRUE,
                  pval.method=TRUE,conf.int = FALSE,
                  font.main=c(20,"bold"),font.x=c(16,"bold"),font.y=c(16,"bold"),
                  font.tickslab=c(14,"plain"), title = tit_female, legend.labs = c("Low", "High"))

splots <- list()
splots[[1]] <- p1m
splots[[2]] <- p1f
p1 <- arrange_ggsurvplots(splots, print = FALSE,
                          ncol = 2, nrow = 1, title = str_c(args,"_",colnames(data)[j]))
png(tit_plot, width=800, height=400)
p1
dev.off()





# ##### To Loop ###
# for(plot_i in 1:length(genes)){
#   list_genes <- genes[plot_i]
#   print(list_genes)
#   ############### Plot variables change genes #########
#   j = grep(list_genes[1],colnames(data))
# 
#   male_sub[,j][male_sub[,j]< median(male_sub[,j])] <- 0
#   male_sub[,j][male_sub[,j]>= median(male_sub[,j])] <- 1
#   female_sub[,j][female_sub[,j]< median(female_sub[,j])] <- 0
#   female_sub[,j][female_sub[,j]>= median(female_sub[,j])] <- 1
# 
#   fit_m <- survfit( Surv(time = PFI.time,
#                          event = PFI)~ male_sub[,j],
#                     data = male_sub)
#   fit_f <- survfit( Surv(time = PFI.time,
#                          event = PFI)~ female_sub[,j],
#                     data = female_sub)
#   tit_male <- str_c(args[1],"_Male_", idmap(colnames(male_sub)[j]))
#   tit_female <- str_c(args[1],"_Female_", idmap(colnames(female_sub)[j]))
#   tit_plot <- str_c(dir, args[1], colnames(data)[j],".png")
#   p1m <- ggsurvplot(fit_m, data=male_sub,palette = c("#404040", "#ca0020"),
#                     ggtheme = theme_bw(),pval = TRUE,
#                     pval.method=TRUE,conf.int = FALSE,
#                     font.main=c(20,"bold"),font.x=c(16,"bold"),font.y=c(16,"bold"),
#                     font.tickslab=c(14,"plain"), title = tit_male, legend.labs = c("Low", "High"))
#   p1f <- ggsurvplot(fit_f, data=female_sub,palette = c("#404040", "#ca0020"),
#                     ggtheme = theme_bw(),pval = TRUE,
#                     pval.method=TRUE,conf.int = FALSE,
#                     font.main=c(20,"bold"),font.x=c(16,"bold"),font.y=c(16,"bold"),
#                     font.tickslab=c(14,"plain"), title = tit_female, legend.labs = c("Low", "High"))
# 
#   splots <- list()
#   splots[[1]] <- p1m
#   splots[[2]] <- p1f
#   p1 <- arrange_ggsurvplots(splots, print = FALSE,
#                             ncol = 2, nrow = 1, title = str_c(args,"_",colnames(data)[j]))
#   png(tit_plot, width=800, height=400)
#   print(p1)
#   dev.off()
# }
