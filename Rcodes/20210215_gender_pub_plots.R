# Barchart and Overlap Barplot_2A.jpeg  Barplot_2B.jpeg  Barplot_S1.jpeg  Barplot_S2.jpeg

##Barchart for number of significant genes with p<0.20 FDR cutoff
setwd("/home/workstation/Documents/Sarthak/gender_prognosis")
#for intfile
library(plyr)
library(readr)
library(stringr)
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore/")
myfiles = list.files(pattern="*_FDR_int.tab", full.names=FALSE)
int_df <- data.frame()
for (i in myfiles){
  int <- read.table(file = i,  sep = '\t', header = TRUE)
  tot <- nrow(int)
  print <- strsplit(i,split="_")[[1]][1]
  int_df <- rbind(int_df, data.frame("Cancer"=print, "Int"=tot))
}

#for estfile
myfiles = list.files(pattern="*_FDR_est.tab", full.names=FALSE)
est_df <- data.frame()
for (i in myfiles){
  est <- read.table(file = i,  sep = '\t', header = TRUE)
  tot <- nrow(est)
  print <- strsplit(i,split="_")[[1]][1]
  est_df <- rbind(est_df, data.frame("Cancer"=print, "Est"=tot))
}
rm(list = setdiff(ls(),c("est_df","int_df")))

#Reading data and initialising it for Barchart for number of significant genes with p<0.20 FDR cutoff 
summary <- merge(est_df, int_df, by = "Cancer")
library(reshape2)
summary.long<-melt(summary,id.vars="Cancer")
library(ggplot2)
p <- ggplot(summary.long,aes(x=Cancer,y=value,fill=factor(variable)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Model",
                      breaks=c("Est", "Int"),
                      labels=c("Without", "With"))+
  xlab("Cancer")+ylab("No of genes")+ ggtitle("Barplot for Cancer vs No of Genes")
p
##Figure 1 ## 
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210215_FDR_gender_zscore")
#jpeg(filename = "Barplot_2A.jpeg")
tiff(filename = "Barplot_2A.tiff", height =1200, width = 1200)
p
dev.off()
rm(list = ls())

######## Calculating overlaps between Int and Est ########
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore")
myfiles1 = list.files(pattern="*_FDR_int.tab", full.names=FALSE)
myfiles2 = list.files(pattern="*_FDR_est.tab", full.names=FALSE)
venn_df <- data.frame()
for(i in 1:length(myfiles1)){
  int <- read.table(file = myfiles1[i],  sep = '\t', header = TRUE)
  int <- subset(int, select = 'gene')
  int$gene <- as.character(int$gene)
  est <- read.table(file = myfiles2[i],  sep = '\t', header = TRUE)
  est <- subset(est, select = 'gene')
  est$gene <- as.character(est$gene)
  Cancer <- strsplit(myfiles1[i],split="_")[[1]][1]
  print(Cancer)
  print(strsplit(myfiles2[i],split="_")[[1]][1])
  Int <- length(int$gene)
  Est <- length(est$gene)
  Overlap <- length(intersect(est$gene, int$gene))
  v <- cbind(Cancer,Int,Est,Overlap)
  venn_df = rbind(venn_df,v)
  print("#####")
}
rm(list = setdiff(ls(),"venn_df"))
venn = venn_df
venn$Int = as.numeric(as.character(venn$Int))
venn$Est = as.integer(as.character(venn$Est))
venn$Overlap = as.integer(as.character(venn$Overlap))
venn$union <- venn$Int+venn$Est-venn$Overlap
venn[,2:5]<- venn[,2:5]/venn[,5]
summary <- venn[,c("Cancer", "Est", "Int", "Overlap")]
library(reshape2)
summary.long<-melt(summary,id.vars="Cancer")
library(ggplot2)
p <- ggplot(summary.long,aes(x=Cancer,y=value,fill=factor(variable)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Model",
                      breaks=c("Est", "Int","Overlap"),
                      labels=c("Without", "With","Overlap"))+
  xlab("Cancer")+ylab("No of genes")+ ggtitle("Barplot for Cancer: Overlap Between Models")
p
##Figure 2 ## 
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210215_FDR_gender_zscore")
#jpeg(filename = "Barplot_S2D.jpeg")
tiff(filename = "Barplot_S2D.tiff", height=1200, width =1200)
p
dev.off()
rm(list = ls())

################ Coefficient Block ######
##These plots were plotted by changing the low/high cutoff to 25,15,10,5
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore")
myfiles = list.files(pattern="*_FDR_int.tab", full.names=FALSE)
int_df <- data.frame()
for (i in myfiles){
  int <- read.table(file = i,  sep = '\t', header = TRUE)
  tot = nrow(int)
  print <- strsplit(i,split="_")[[1]][1]
  int_low <- nrow(int[int$Interaction_exp.coef.<0.75,])
  int_high <- nrow(int[int$Interaction_exp.coef.>1.25,])
  int_df <- rbind(int_df, data.frame("Cancer"=print, "Significant"=tot,
                                     "Low_Coeff"=int_low, "High_Coeff"=int_high))
}
rm(list = setdiff(ls(),"int_df"))
max_range = ceiling(max(int_df[,c(3,4)])/100)*100
summary <- int_df
p2a <-ggplot(summary, aes(Cancer, Low_Coeff)) 
p2a <- p2a + geom_bar(stat = "identity", position = "dodge", fill = "#0073C2FF") 
p2a <- p2a + xlab("CancerType") + ylab("No. of Genes") + ggtitle("Genes decreasing the risk")
p2a <- p2a +ylim(c(0,max_range)) + coord_flip() 
p2a <- p2a + scale_y_reverse()
# high coefficient 
p2b <-ggplot(summary, aes(Cancer, High_Coeff)) 
p2b <- p2b + geom_bar(stat = "identity", position = "dodge", fill = "#0073C2FF") 
p2b <- p2b + xlab("CancerType") + ylab("No. of Genes") + ggtitle("Genes increasing the risk")
p2b <- p2b +ylim(c(0,max_range))+coord_flip()
library(gridExtra)

grid.arrange(
  p2a,
  p2b,
  nrow = 1,
  top = "Genes that affect survival >= 25%"
)
##Figure 3 ## 
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210215_FDR_gender_zscore")
#jpeg(filename = "Barplot_2B.jpeg")
#jpeg(filename = "Barplot_S3B.jpeg")
tiff(filename = "Barplot_S3B.tiff", width = 1200, height =1200)
grid.arrange(
  p2a,
  p2b,
  nrow = 1,
  top = "Genes that affect survival >= 50%"
)

dev.off()

##### Common genes #####
rm(list= ls())
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore_name")
count_p <- read.csv("new.txt", sep="")
library(dplyr)
df <- count_p %>%
  group_by(frequency) %>%
  summarise(counts = n())
library(ggplot2)
p <-ggplot(df, aes(x = frequency, y = counts)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) + xlab("Number of cancer types") +
  ggtitle("Distribution of significant genes across Cancer types")
p
##Figure 4 ## 
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/plots/20210215_FDR_gender_zscore")
#jpeg(filename = "Barplot_2B.jpeg")
tiff(filename = "Barplot_2B.tiff", width = 1200, height = 1200)
p
dev.off()
