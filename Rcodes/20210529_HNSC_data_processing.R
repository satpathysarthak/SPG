## Microarray data processing
library(data.table)
setwd('/home/workstation/Documents/Sarthak/gender_prognosis/validation/HNSC')
# sample <- read.delim("~/Documents/Sarthak/gender_prognosis/validation/HNSC/sample.txt", header=FALSE, stringsAsFactors = F)
# sample <- data.table(t(sample))
# colnames(sample) <- as.character(sample[1,])
# sample <- sample[-1,]
# write.table(sample, 'output/sample.tab', append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
# tum <- fread('GSE65858_Non-normalized_data.txt', header = T)
# tum <- data.frame(tum)
# tum <- tum[-grep('Detection.Pval',colnames(tum))]
# rownames(tum) <- tum$ID_REF
# tum$ID_REF <- NULL
# t_row <- rownames(tum)
# t_col <- colnames(tum)
# tum <- transpose(tum)
# rownames(tum) <- t_col
# colnames(tum) <- t_row
# normalize <- function(df, cols) {
#   result <- df # make a copy of the input data frame
#   
#   for (j in cols) { # each specified col
#     m <- mean(df[,j]) # column mean
#     std <- sd(df[,j]) # column (sample) sd
#     
#     for (i in 1:nrow(result)) { # each row of cur col
#       result[i,j] <- (result[i,j] - m) / std
#     }
#   }
#   return(result)
# }
# tum <- normalize(tum,colnames(tum))
# write.table(tum, 'output/tum_norm.tab', append = FALSE, col.names =TRUE,row.names = T, quote = T, sep='\t')
sample = fread('output/sample.tab', header=T)
tum <- fread('output/tum_norm.tab',header = T)
tum <- data.frame(tum)
rownames(tum) <- tum$ID
tum$ID <- NULL
t_row <- rownames(tum)
t_col <- colnames(tum)
tum <- transpose(tum)
rownames(tum) <- t_col
colnames(tum) <- t_row
tum$probe <- rownames(tum)
#df_name <- data.frame('probe'=rownames(tum))
map <- read.table('output/gene_map.txt', header = TRUE, sep = "\t", 
           col.names = paste0("V",seq_len(2)), fill = TRUE, stringsAsFactors = F)
colnames(map) <- c('probe','gene')
map$gene[map$gene==""]<-NA
map <- unique.data.frame(map)
tum_mapped <- merge(tum,map, by='probe', all.x=T, all.y = F)
tum_mapped <- data.table(tum_mapped)
tum_mapped[, `:=` (count = .N), by = gene]
tum_mapped <- data.frame(tum_mapped)
tum_mapped$median <- NA
for(i in 1:nrow(tum_mapped)){
  tum_mapped$median[i] <- median(as.numeric(tum_mapped[i,colnames(tum_mapped)[grep('GW',colnames(tum_mapped))]]))
}
tum_mapped <- data.table(tum_mapped)
tum_mapped[, `:=` (maxmed = max(median)), by = gene]
tum_final <- tum_mapped[median-maxmed==0]
tum_final$count = NULL
tum_final$median = NULL
tum_final$probe = NULL
tum_final$maxmed = NULL
tum <- data.frame(tum_final)
tum_mapped = NULL
tum_final = NULL
tum <- na.omit(tum)
rownames(tum) <- tum$gene
tum$gene = NULL
t_row <- rownames(tum)
t_col <- colnames(tum)
tum <- transpose(tum)
rownames(tum) <- t_col
colnames(tum) <- t_row
tum$ID <- rownames(tum)

sample$gender <- as.character(sample$gender)
sample[sample$gender=="F",]$gender <- "1"
sample[sample$gender=="M",]$gender <- "0"
sample$gender <- as.numeric(sample$gender)

data <- merge(sample, tum, by="ID", all = F)
data <- data.frame(data)
write.table(data, 'output/data_HNSC.tab', append = FALSE, col.names =TRUE,row.names = FALSE, quote = FALSE, sep='\t')
