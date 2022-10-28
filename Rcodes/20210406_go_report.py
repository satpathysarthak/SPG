import pandas as pd
import sys
indir = "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/reports/"
cancer = sys.argv[1]
#cancer='BLCA_ns_s'
infile = indir + cancer + '_report.xlsx'
xls = pd.ExcelFile(infile)
df1 = pd.read_excel(xls, 'Molecular function', skiprows=8, header =0)
df2 = pd.read_excel(xls, 'Biological pathway', skiprows=8, header =0)
df3 = pd.read_excel(xls, 'Transcription factor', skiprows=8, header =0)
df4 = pd.read_excel(xls, 'Clinical phenotype', skiprows=8, header =0)
df1['cancer_sig']=cancer
df2['cancer_sig']=cancer
df3['cancer_sig']=cancer
df4['cancer_sig']=cancer
outdir1 = "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/reports_all/bhfilter_"
out1_df1 = df1[df1['BH method']<0.1]
out1_df2 = df2[df2['BH method']<0.1]
out1_df3 = df3[df3['BH method']<0.1]
out1_df4 = df4[df4['BH method']<0.1]
out1_df1.to_csv((outdir1+'MF.txt'),index=False,header = False,sep="\t", mode ='a')
out1_df2.to_csv((outdir1+'BP.txt'),index=False,header = False,sep="\t",mode ='a')
out1_df3.to_csv((outdir1+'TF.txt'),index=False,header = False,sep="\t", mode='a')
out1_df4.to_csv((outdir1+'CP.txt'),index=False,header = False,sep="\t",mode = 'a')
outdir2 = "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/reports_all/top10_"
out2_df1 = df1[df1['P-value (Hypergeometric test)']<0.05]
out2_df1 = out2_df1.nlargest(10, ['No. of genes in the dataset'])
out2_df2 = df2[df2['P-value (Hypergeometric test)']<0.05]
out2_df2 = out2_df2.nlargest(10, ['No. of genes in the dataset'])
out2_df3 = df3[df3['P-value (Hypergeometric test)']<0.05]
out2_df3 = out2_df3.nlargest(10, ['No. of genes in the dataset'])
out2_df4 = df4[df4['P-value (Hypergeometric test)']<0.05]
out2_df4 = out2_df4.nlargest(10, ['No. of genes in the dataset'])
out2_df1.to_csv((outdir2+'MF.txt'),index=False,header = False,sep="\t", mode ='a')
out2_df2.to_csv((outdir2+'BP.txt'),index=False,header = False,sep="\t",mode ='a')
out2_df3.to_csv((outdir2+'TF.txt'),index=False,header = False,sep="\t", mode='a')
out2_df4.to_csv((outdir2+'CP.txt'),index=False,header = False,sep="\t",mode = 'a')