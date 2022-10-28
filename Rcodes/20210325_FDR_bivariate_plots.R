library(data.table)
library(ggplot2)
library(ggrepel)
input_dir <- "/home/sarthak/Documents/projects/MKlab/cox_regression/bivariate_gender/"
input_files_male  <- list.files(input_dir,pattern = "_male.tab",full.names = F)
input_files_female  <- list.files(input_dir,pattern = "_female.tab",full.names = F)
outdir <- "/home/sarthak/Documents/projects/MKlab/plots/20210324_bivariate/"
setwd(input_dir)
map_dir = "/home/sarthak/Documents/projects/MKlab/cox_regression/FDR_gender_zscore/withnames/"
input_files_map <- list.files(map_dir,pattern = "_map_int.txt",full.names = T)
data_old <- data.frame()
for(file_i in 1:length(input_files_male)){
  map <- fread(input_files_map[file_i], header = F)
  colnames(map) <- c("gene","gene_name")
  data_male <- fread(input_files_male[file_i])
  data_male <- data_male[,c(1:6)]
  data_female <- fread(input_files_female[file_i])
  data_female <- data_female[,c(1:6)]
  colnames(data_male)[c(2:6)] <- gsub("gene","male_gene",colnames(data_male)[c(2:6)])
  colnames(data_female)[c(2:6)] <- gsub("gene","female_gene",colnames(data_female)[c(2:6)])
  data <- merge(data_male,data_female,by='gene')
  data$direction_male <- ifelse(data$male_gene_Pr...z..<0.05,1,0)
  data$direction_female <- ifelse(data$female_gene_Pr...z..<0.05,1,0)
  print(table(data$direction_male,data$direction_female))
  cancer <- sapply(strsplit(input_files_male[file_i], split = "_",fixed = TRUE),"[",1)
  dir <- as.vector(table(data$direction_male,data$direction_female))
  data$sig <- NA
  data$sig <- as.character(paste0(data$direction_male,data$direction_female))
  data$sig_dir <- data$male_gene_coef*data$female_gene_coef
  data$sig_dir <- data$sig_dir/abs(data$sig_dir)
  data[data$sig == "00",]$sig <- "ns_ns"
  data[data$sig == "01",]$sig <- "ns_s"
  data[data$sig == "10",]$sig <- "s_ns"
  data[data$sig == "11" & data$ sig_dir == "-1",]$sig <- "s_s_opp"
  data[data$sig == "11" & data$ sig_dir == "1",]$sig <- "s_s_same"
  data <- merge(data,map, by='gene')
  
  per_male <- quantile(abs(data$male_gene_coef),0.50)
  per_female <- quantile(abs(data$female_gene_coef),0.50)
  colorset = c('ns_ns'='#404040','ns_s'='#f77ec8','s_ns'='#6fbde3','s_s_opp'='red','s_s_same'='darkgrey')
  sp <- ggplot(data = data, aes(x=male_gene_coef,y=female_gene_coef,color=sig))+
    geom_point()+ xlab(paste0("Male Beta")) + ylab("Female Beta") + ggtitle(cancer)
  # sp <- sp+scale_color_manual(values=c("#404040","#F1B6DA","#92C5DE","red", "darkgrey"))
  sp <- sp+scale_color_manual(values=colorset)
  sp <- sp+theme_classic(base_size = 22)+geom_hline(yintercept=0,linetype = "dotted")+geom_vline(xintercept=0,linetype = "dotted")
  #sp <- sp + xlim(c(lab_all*-1,lab_all)) + ylim(c(lab_all*-1,lab_all))
  sp <- sp + guides(color = guide_legend(override.aes = list(size = 2)))
  #colorset = c('ns_ns'='#404040','ns_s'='#F1B6DA','s_ns'='#92C5DE','s_s_opp'='green','s_s_same'='darkgrey')
  sp <- sp + geom_text_repel(aes(label=ifelse(sig%in%c('ns_s','s_ns','s_s_opp')&(abs(male_gene_coef)>per_male |abs(female_gene_coef)>per_female),as.character(gene_name),'')),size=4)
  tiff(paste0(outdir,cancer,"_beta.tiff"),height=1200,width=1834,res=200)
  print(sp)
  dev.off()
}
