## KM plot specific genes
list_genes = c("ENSG00000155850.7","ENSG00000172493.19","ENSG00000158711.12")
library(ggplot2)
args<-"HNSC"
#### Input & Output file names ####
tumor_file = str_c("/media/sarthak/gender_prognosis/Tumor/",args[1],".tab")
clinical_file <- str_c("/media/sarthak/gender_prognosis/Clinical/processed/",args[1],"_pro.tab")
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

male_sub <- data[data$gender==0,]
female_sub <- data[data$gender==1,]
data[data$gender==1,]$gender <- "female"
data[data$gender==0,]$gender <- "male"
j = grep(list_genes[3],colnames(data))
name <- colnames(data)[j]
name <- strsplit(name, split = ".",fixed = T)
name = name[[1]][1]
p<-ggplot(data, mapping = aes_string(x = names(data)[6], y = names(data)[j])) +
  geom_boxplot()
p
tit_plot <- str_c("/media/sarthak/gender_prognosis/plots/", name,"_barplot_exp.png")
png(tit_plot, width=800, height=400)
print(p)
dev.off()
