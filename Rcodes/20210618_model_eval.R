## baplots
library(data.table)
library(ggplot2)
setwd('/home/workstation/Documents/Sarthak/gender_prognosis/')
data <- fread('analysis/prognostic/HNSC_AUC.tab', header = T)
data$TCGA_diff <- data$TCGA_AUC_both_all - data$TCGA_AUC_est_all
data$GEO_diff <- data$GEO_AUC_both_all - data$GEO_AUC_est_all
df_long <- subset(data, select = c(1,14,15))
df_long <- melt(df_long, id.vars=c('bootstrap'))
df_long$variable <- gsub('_diff','',df_long$variable)
ggplot(data = df_long, aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))

df_long <- subset(data, select = c(1,14,15))
p <- ggplot(data, aes(x = reorder(bootstrap, -TCGA_diff), y = TCGA_diff))
p <- p + geom_bar(stat="identity", color='skyblue',fill='steelblue')
p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
p
p <- ggplot(data, aes(x = reorder(bootstrap, -GEO_diff), y = GEO_diff))
p <- p + geom_bar(stat="identity", color='skyblue',fill='steelblue')
p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
p

df_long <- subset(data, select = c(1,14,15))
df_long <- df_long[order(-df_long$TCGA_diff),]
df_long$col_ord <- rownames(df_long)
df_long <- melt(df_long, id.vars=c('bootstrap','col_ord'))
df_long$variable <- gsub('_diff','',df_long$variable)
df_long$bootstrap <- as.factor(df_long$bootstrap)
df_long$variable <- factor(df_long$variable, levels =c('TCGA','GEO'))
p <- ggplot(df_long, aes(x = reorder(bootstrap, col_ord),y = value, fill = variable))
p <- p + geom_bar(stat="identity", position=position_dodge())
p
