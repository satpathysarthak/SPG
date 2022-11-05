# Plots from Enrichment scores
library(topGO)
library(stringr)
#library(data.table)
output_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/plots/"
dir_path <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/"
dir <- list.files(dir_path, pattern= "_same.rnk", all.files=FALSE,  full.names=FALSE)
setwd(dir_path)
for(file_name in dir){
#for(file_name in c("KIRP_s_s_opp.rnk","THCA_s_s_opp.rnk")){
  print(file_name)
  id_file <- file_name
  file_name <- strsplit(file_name,".",fixed = T)[[1]][1]
  id <- read.delim(id_file, header=FALSE, comment.char="#")
  if(nrow(id)<20){next}
  id$V2 <- 0.01
  genes <- id$V2
  names(genes) <- id$V1
  selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
  #id$V2 <- NULL
  allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
  GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=10)
  results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
  #if(length(results.ks)==0){next}
  goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
  goEnrichment$KS <- gsub("< ","",goEnrichment$KS)
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
  goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
  goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
  require(ggplot2)
  goEnrichment <- na.omit(goEnrichment)
  p <- ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("Biological Processes") +
    ylab("Enrichment") +
    ggtitle(file_name) +
    scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
    theme_bw(base_size=24) +
    theme(
      legend.position='none',
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
      axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
      axis.title=element_text(size=24, face="bold"),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=18),  #Text size
      title=element_text(size=18)) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
  jpeg(str_c(output_dir,file_name,".jpg"),width = 930,height = 539)
  print(p)
  dev.off()
  #rm(list = c("results.ks","goEnrichment","p","GOdata"))
}

