## Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated
rm(list= ls())
#args<-commandArgs(TRUE)
args = 'GBM'
print(args)
setwd("/home/sarthak/Documents/projects/MKlab/")
library(data.table)
library(stringr)
#### Input & Output file names ####
tumor_file = str_c("/home/sarthak/Documents/projects/MKlab/Tumor/",args[1],".tab")
clinical_file <- str_c("/home/sarthak/Documents/projects/MKlab/Clinical/processed/",args[1],"_pro.tab")
#fdr_file <- str_c("/home/sarthak/Documents/projects/MKlab/cox_regression/FDR_gender_zscore_name/",args[1],"_FDR_int.name")
fdr_file <- str_c("/home/sarthak/Documents/projects/MKlab/analysis/bivariate_go/","GBM_s_s_opp.tab")
out_file <- str_c("/home/sarthak/Documents/projects/MKlab/plots/20210310_heatmaps/",args[1],"_heatmap_IDH.pdf")
#output_file_est <- str_c("/home/sarthak/Documents/projects/MKlab/cox_regression/gender_zscore/",args[1],"_est.tab")
#output_file_int <- str_c("/home/sarthak/Documents/projects/MKlab/cox_regression/gender_zscore/",args[1],"_int.tab")
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
#sample$submitter_id = sample$V1
sample$ID <- sapply(strsplit(sample$submitter_id,split = "-"), "[",3)
sample$gender <- as.character(sample$gender)
sample[sample$gender=="female",]$gender <- "1"
sample[sample$gender=="male",]$gender <- "0"
sample$gender <- as.numeric(sample$gender)
sample <- sample[!is.na(sample$gender),]
### sample2
clinical_file <- str_c("/home/sarthak/Documents/projects/MKlab/Clinical/Clinical_OS/",args[1],".tab")
sample2 = fread(clinical_file)
sample2 = sample2[,c(1,6)]
foo = merge(x=sample,y=sample2,by.x="submitter_id",by.y = "V1", all.x=T, ally.y=F)
foo$IDH_status_WT = foo$`IDH status_WT`
foo$`IDH status_WT` = NULL
foo = na.omit(foo)
sample = foo
##### Data Merging #####
data <- merge(sample, tum, by="ID", all = F)
data <- data.frame(data)
##### Survival analysis: Cox proportional Hazard Model
#data <- data[order(data$gender),]
data <- data[order(data$IDH_status_WT),]
### BLOCK A: Simple heatmaps
# BiocManager::install("ComplexHeatmap")
data[data$IDH_status_WT==1,]$IDH_status_WT<- "Mut"
data[data$IDH_status_WT==0,]$IDH_status_WT <- "WT"
library(ComplexHeatmap)
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,8:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames
col = list(IDH_status_WT = c("Mut"="#40B0A6","WT"="#E1BE6A"))
# """
# clinical_file <- str_c("/home/sarthak/Documents/projects/MKlab/Clinical/Clinical_OS/",args[1],".tab")
# sample <- fread(clinical_file)
# foo = merge(x=data,y=sample,by.x="submitter_id",by.y = "V1", all.x=T, ally.y=F)
# """

# Create the heatmap annotation
out_file <- str_c("/home/sarthak/Documents/projects/MKlab/plots/20210310_heatmaps/",args[1],"_heatmap_IDH_s_s_opp.pdf")
ha <- HeatmapAnnotation(
  IDH_status_WT = data$IDH_status_WT,
  col = col
)
p <- Heatmap(t(mat_data),
             name = "z-score", #title of legend
             column_title = str_c(args[1]," Samples") , row_title = "Genes",
             #row_names_gp = gpar(fontsize = 7), # Text size for row names
             top_annotation = ha, cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F
)
  pdf(out_file)
  print(p)
  dev.off()

#### BLOCK B: heatmaps with median values

# data[data$gender=="1",]$gender <- "female"
# data[data$gender=="0",]$gender <- "male"
# data <- data[,-c(1:5)]
# data$gender <- as.character(data$gender)
# gender = data$gender
# data <-data.table(data)
# data <- data[,lapply(.SD,median),by=gender]
# data <- data.frame(data)
# rnames <- data[,1]                            # assign labels in column 1 to "rnames"
# mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
# rownames(mat_data) <- rnames
# col = list(gender = c("male"="blue","female"="green"))
# # Create the heatmap annotation
# library(ComplexHeatmap)
# ha <- HeatmapAnnotation(
#   gender = data$gender,
#   col = col
# )
# p <- Heatmap(t(mat_data),
#         name = "z-score", #title of legend
#         column_title = str_c(args[1]," Samples") , row_title = "Genes",
#         row_names_gp = gpar(fontsize = 7), # Text size for row names
#         #top_annotation = ha, cluster_rows = F,cluster_columns = F,show_column_names = F
#         cluster_rows = F,cluster_columns = F,show_column_names = T
# )
# jpeg(out_file,)
# print(p)
# dev.off()

