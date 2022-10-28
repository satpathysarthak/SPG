##Barchart for number of significant genes with p<0.20 FDR cutoff
setwd("/home/sarthak/Documents/projects/MKlab/")
#for intfile
library(plyr)
library(readr)
library(stringr)
library(ggplot2)
setwd("/home/sarthak/Documents/projects/MKlab/analysis/prelim/")
myfiles = list.files(pattern="*.tab", full.names=FALSE)
outdir <- '/home/sarthak/Documents/projects/MKlab/plots/20210618_prelim/'
int_df <- data.frame()
for (i in myfiles){
  int <- read.table(file = i,  sep = '\t', header = TRUE)
  tot <- nrow(int)
  print <- strsplit(i,split="_")[[1]][1]
  m = table(int$gender)[1]
  f = table(int$gender)[2]
  cen = sum(int$PFI==0)
  uncen = sum(int$PFI==1)
  #int_df <- rbind(int_df, data.frame("Cancer"=print, "total"=tot, 'm'=m, 'f'=f))
  int_df <- rbind(int_df, data.frame("Cancer"=print, "total"=tot, 'cen'=cen, 'uncen'=uncen))
}


## Barcharts with ceeeeeeeeensored uncensored information
int_df <- data.frame()
for (i in myfiles){
  int <- read.table(file = i,  sep = '\t', header = TRUE)
  print <- strsplit(i,split="_")[[1]][1]
  int$cancer <- print
  int <-   subset(int, select= c("cancer", "submitter_id", "PFI", "age", "gender"))
  int_df <- rbind(int_df, int)
}
data <- data.frame(int_df)
df_long <- reshape2::melt(data, id.vars=c("cancer","submitter_id","age"))
df_long = data
df_long$gender <- ifelse(df_long$gender == 0, "M","F")
df_long$gender <- factor(df_long$gender, levels = c("M","F"))
df_long$PFI <- as.factor(df_long$PFI)
p <- ggplot(data=df_long, aes(x= gender,fill=PFI))
p <- p + geom_bar(stat = "count")+facet_grid(~cancer)
p
outfile = str_c(outdir,"gender.tiff")
tiff(outfile,height=685,width=1200)
print(p)
dev.off()

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
          "#F4EDCA", "#52854C", "#4E84C4", "#293352")

p <- ggplot(data, aes(x=age, colour=cancer)) + geom_density() + scale_fill_manual(cbp2)
outfile = str_c(outdir,"age.tiff")
tiff(outfile,height=685,width=1200)
print(p)
dev.off()

