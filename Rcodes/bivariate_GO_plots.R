# Plots from processed GO analysis
dir_path <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/reports_all'
file_list <- list.files(dir_path, pattern= "BP.sub", all.files=FALSE,  full.names=FALSE)
setwd(dir_path)
i=2
data <- fread(file_list[i])
data$FDR_BH_p = -log10(data$FDR_BH_p)
name = strsplit(file_list[i], ".", fixed =T)[[1]][1]
ns_ns = subset(data, sig == 'ns_ns')
ns_s = subset(data, sig == 'ns_s')
s_ns = subset(data, sig == 's_ns')
s_s_opp = subset(data, sig == 's_s_opp')
s_s_same = subset(data, sig == 's_s_same')
# ns_ns
ns_ns = ns_ns[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Biological_pathway]
p <- ns_ns %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Biological_pathway,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': ns_ns')) +
  xlab("avg -log10(pval)") + ylab("Biological Pathway")
tiff(paste0(name, '_ns_ns.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()
# ns_s
ns_s = ns_s[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Biological_pathway]
p <- ns_s %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Biological_pathway,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': ns_s')) +
  xlab("avg -log10(pval)") + ylab("Biological Pathway")
tiff(paste0(name, '_ns_s.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_ns
s_ns = s_ns[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Biological_pathway]
p <- s_ns %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Biological_pathway,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_ns')) +
  xlab("avg -log10(pval)") + ylab("Biological Pathway")
tiff(paste0(name, '_s_ns.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_s_opp
s_s_opp = s_s_opp[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Biological_pathway]
p <- s_s_opp %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Biological_pathway,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_s_opp')) +
  xlab("avg -log10(pval)") + ylab("Biological Pathway")
tiff(paste0(name, '_s_s_opp.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_s_same
s_s_same = s_s_same[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Biological_pathway]
p <- s_s_same %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Biological_pathway,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_s_same')) +
  xlab("avg -log10(pval)") + ylab("Biological Pathway")
tiff(paste0(name, '_s_s_same.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()


##################################

dir_path <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/reports_all'
file_list <- list.files(dir_path, pattern= "MF.sub", all.files=FALSE,  full.names=FALSE)
setwd(dir_path)
i=2
data <- fread(file_list[i])
name = strsplit(file_list[i], ".", fixed =T)[[1]][1]
data$FDR_BH_p = -log10(data$FDR_BH_p)
ns_ns = subset(data, sig == 'ns_ns')
ns_s = subset(data, sig == 'ns_s')
s_ns = subset(data, sig == 's_ns')
s_s_opp = subset(data, sig == 's_s_opp')
s_s_same = subset(data, sig == 's_s_same')
# ns_ns
ns_ns = ns_ns[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Molecular_function]
p <- ns_ns %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Molecular_function,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': ns_ns')) +
  xlab("avg -log10(pval)") + ylab("Molecular function")
tiff(paste0(name, '_ns_ns.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()
# ns_s
ns_s = ns_s[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Molecular_function]
p <- ns_s %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Molecular_function,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': ns_s')) +
  xlab("avg -log10(pval)") + ylab("Molecular function")
tiff(paste0(name, '_ns_s.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_ns
s_ns = s_ns[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Molecular_function]
p <- s_ns %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Molecular_function,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_ns')) +
  xlab("avg -log10(pval)") + ylab("Molecular function")
tiff(paste0(name, '_s_ns.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_s_opp
s_s_opp = s_s_opp[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Molecular_function]
p <- s_s_opp %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Molecular_function,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_s_opp')) +
  xlab("avg -log10(pval)") + ylab("Molecular function")
tiff(paste0(name, '_s_s_opp.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_s_same
s_s_same = s_s_same[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Molecular_function]
p <- s_s_same %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Molecular_function,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_s_same')) +
  xlab("avg -log10(pval)") + ylab("Molecular function")
tiff(paste0(name, '_s_s_same.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()


##############################

dir_path <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/reports_all'
file_list <- list.files(dir_path, pattern= "TF.sub", all.files=FALSE,  full.names=FALSE)
setwd(dir_path)
i=2
data <- fread(file_list[i])
name = strsplit(file_list[i], ".", fixed =T)[[1]][1]
data$FDR_BH_p = -log10(data$FDR_BH_p)
ns_ns = subset(data, sig == 'ns_ns')
ns_s = subset(data, sig == 'ns_s')
s_ns = subset(data, sig == 's_ns')
s_s_opp = subset(data, sig == 's_s_opp')
s_s_same = subset(data, sig == 's_s_same')
# ns_ns
ns_ns = ns_ns[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Transcription_factor]
p <- ns_ns %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Transcription_factor,avg_p),fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': ns_ns')) +
  xlab("avg -log10(pval)") + ylab("Transcription factor")
tiff(paste0(name, '_ns_ns.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()
# ns_s
ns_s = ns_s[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Transcription_factor]
p <- ns_s %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Transcription_factor,avg_p),fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': ns_s')) +
  xlab("avg -log10(pval)") + ylab("Transcription factor")
tiff(paste0(name, '_ns_s.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_ns
s_ns = s_ns[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Transcription_factor]
p <- s_ns %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Transcription_factor,avg_p),fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_ns')) +
  xlab("avg -log10(pval)") + ylab("Transcription factor")
tiff(paste0(name, '_s_ns.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_s_opp
s_s_opp = s_s_opp[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Transcription_factor]
p <- s_s_opp %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Transcription_factor,avg_p),fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_s_opp')) +
  xlab("avg -log10(pval)") + ylab("Transcription factor")
tiff(paste0(name, '_s_s_opp.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_s_same
s_s_same = s_s_same[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Transcription_factor]
p <- s_s_same %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Transcription_factor,avg_p),fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_s_same')) +
  xlab("avg -log10(pval)") + ylab("Transcription factor")
tiff(paste0(name, '_s_s_same.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()


#############

dir_path <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/reports_all'
file_list <- list.files(dir_path, pattern= "CP.sub", all.files=FALSE,  full.names=FALSE)
setwd(dir_path)
i=2
data <- fread(file_list[i])
name = strsplit(file_list[i], ".", fixed =T)[[1]][1]
data$FDR_BH_p = -log10(data$FDR_BH_p)
ns_ns = subset(data, sig == 'ns_ns')
ns_s = subset(data, sig == 'ns_s')
s_ns = subset(data, sig == 's_ns')
s_s_opp = subset(data, sig == 's_s_opp')
s_s_same = subset(data, sig == 's_s_same')
# ns_ns
ns_ns = ns_ns[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Clinical_phenotype]
p <- ns_ns %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Clinical_phenotype,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': ns_ns')) +
  xlab("avg -log10(pval)") + ylab("Clinical phenotype")
tiff(paste0(name, '_ns_ns.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()
# ns_s
ns_s = ns_s[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Clinical_phenotype]
p <- ns_s %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Clinical_phenotype,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': ns_s')) +
  xlab("avg -log10(pval)") + ylab("Clinical phenotype")
tiff(paste0(name, '_ns_s.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_ns
s_ns = s_ns[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Clinical_phenotype]
p <- s_ns %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Clinical_phenotype,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_ns')) +
  xlab("avg -log10(pval)") + ylab("Clinical phenotype")
tiff(paste0(name, '_s_ns.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_s_opp
s_s_opp = s_s_opp[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Clinical_phenotype]
p <- s_s_opp %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Clinical_phenotype,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_s_opp')) +
  xlab("avg -log10(pval)") + ylab("Clinical phenotype")
tiff(paste0(name, '_s_s_opp.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

# s_s_same
s_s_same = s_s_same[, .(count = .N, avg_p = mean(FDR_BH_p)), by = Clinical_phenotype]
p <- s_s_same %>%  
  arrange(desc(avg_p))%>%
  head(30)%>%
  ggplot(aes(x=avg_p, y=reorder(Clinical_phenotype,avg_p), fill=count))+
  geom_bar(stat="identity") + ggtitle(paste0(name, ': s_s_same')) +
  xlab("avg -log10(pval)") + ylab("Clinical phenotype")
tiff(paste0(name, '_s_s_same.tiff'),height=1200,width=1834,res=200)
print(p)
dev.off()

