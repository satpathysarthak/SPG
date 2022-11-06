# SPG

The repository contains Scripts and Processed Data for SPG manuscript.

## Cox Regression model results

cox_regression/ has all the cox regression models.
cox_regression/FDR_gender_zscore/ has the cox-regression model which 

## Scripts
bivariate_GO_plots.R
 Plots from processed GO analysis

cox_gender_miR.R
Script automated for miR gender Cox proportional regression

cox_gender.R
 Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated

cox_gender_z.R
 Script to perform Survival Analysis with Tumor data post z-score transformation using 2 models: Interaction and Estimated

cox_menopausal.R
 Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated

cox_uni_gender.R
 Script to perform Survival Analysis with Tumor data using univariate models on subdata

create_processed_sample.R
 Merge information on gender with Clinical-OS: Creating processed Clinical Sample File

FDR_bivariate_GO_create.R
 Creates a ranked file from output of FunrichR

FDR_bivariate_plots.R
 Plot female vs male beta for ns_ns,s_ns,ns_s,s_s

FDR_bivariate.R
 Script to perform Survival Analysis with Tumor data using 2 models: Male Female

FDR_bivariate_table.R
 Genes in the category: ns_ns,s_ns,ns_s,s_s

FDR_gender_miR.R
 miRNA regression: FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 

FDR_gender.R
 FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.10

FDR_gender_zscore.R
 FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.10

FDR_plots.R
 Plot for miRNA and mRNA expreesion: FDR Filter to get regression results which qualify the mean absolute cutoff and FDR

FDR_plots_zscore.R
 FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.10. Plot the changes. 

FDR_uni_gender.R
 FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.20

gender_pub_plots.R
 Barchart and Overlap Barplot_2A.jpeg  Barplot_2B.jpeg  Barplot_S1.jpeg  Barplot_S2.jpeg

go_plots.R
 Plots from Enrichment scores

go_report.py
 Processing the reports generated by FunrichR

go_rnk_file_creation.R
ID mapping: Specific from taking input from .int file

gsea_rnk_file_creation.R
ID mapping: Specific from taking input from .int file

heatmap_IDH.R
 Heatmap for gene expression betwwen IDH status

heatmap.R
 Heatmap of gene expression for genes from Survial Analysis

hist_sample.R
 Script to generate plots for Gene Expression

HNSC_data_processing.R
 Microarray data processing

HNSC_validation_AUC.R
 Calculates AUC for HNSC validation data

HNSC_validation_KM.R
 KM  plots for HNSC validation and training data

idmap_fun.RData
\8B\00\00\00\00\00\00\B5X\FFk\DBF?\CB\DF\E5\A6\F1\976Q\93lӼ\D0\DA\E0Ą\FC2J\CDH\C9<B\B7P\9A\AC\C20\8E-'Jl\C9H2Y	\85\FD\C9c\BF\8EŻ\D3=\D9\CF'y\B1o\89\E1x_\EE\EE\F3\DE\E9\BD{z\F2\87\C3O\FB\EA'\95'\F1Iʒ䯧͝\EF	I(T\88\91\C92\A5\D9\B4\87\84\D0U\B3)\D7\E98F\8FB\ACP)\E3\CF\C205\C4\C7`h\C2\EE4\DD\DD3\FBe\84\FF\FE-\F7M\CBp)\93\9F*\E3\DF\FD͎\DEYϴ\AD\8A\D5-\C3r\8D\C1y\BFz\AB\EA:V\B0\95n{8\EC\AE\B8\9E\E3\FB\A67\B3\A1V\DE-\D7z\E6\EFF\B7qZ\AD\95\CFʵ\BD*\850\AD\E1\C8k1[􌭞cR\B9~i\8C\FA\8D\ED\\BB^\9BY\AFڝ\D1\C0\B0<\B7~\D2v\BC\CB\F6u\FD°\BA\86\D3:\F6\85e\BB\A6[?:\FC\85z`Zu\DF\F9\C0\F4\AE\D7>/\87my\F6\C3Y\FAɰ\8C\D61=n`\8B`\E8\8E\D1\EE\EEv\8D\BE9\A8D\9C\B4\A6_\D2y\C3i4~>\F9\B1F\9FEt\DCf\BB\E3َ\AB7\F4f5x\CCˤr\F9\F8\80\A1\F7l\A7rEϣｶ\FB\A6B\AB~D\E8!{L\DA\FE\B8v\F5\9B\DEh衠\B1B\F4\EF\D1\C50uN]\BE\F6\F9/*\F7ų\979}\96\C6R'a\CE5\81˳g\00\B0ɉC\87\D0\F5\99\B3PL\B4\8F{Q\BD\C5P\C7\FC\CC\DEȱ\A6k\AB\EAU\B85v~\B6`\F6\E2\F8\BF;>foY\A7\DFv\83[\A6\802\B7\B3c?\87/,@\A4ǂ\89IE\91cOpة\FC7_\A7\A4\F8:\85!\96\D8:\A1\AA\E8\88\EAH\8E\F1}q\BA/1^ī	\AA\E2\D7>N\DF\FFdE5\8EƵ@\E3Ҩ	@H\A8	>$Q\93\80\C4\E8;\E0߁\9C\94FM\A3\E0 \A7\A4QӀ\F4Dx\C6i\AE\93D]F\8F\80?yE\F5) 1Z\BE\F2Si\D4U@z.<\81U\AE\93D]F_\FF

idmap_KM.R
ID mapping: Specific from taking input from .int file

idmap_one.R
idmap <- function(name_ensembl){

idmap_table.R
input_file_ref_from <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_ensembl.tab"

int_beta_scatter.R
 Creates a scatter plot male and female coefficients 

int_est_overlap.R
 Calculating overlap between Main-effect and Interaction model 

KM_genes.R
 KM plot specific genes

KM_plots.R
 Script for plotting KM Plots`

list_gender.R
 Script to generate plots for data exploration: ditribution for cencored and gender data`

ListOfSigNS.R
 Counting ns_ns,s_ns,ns_s,s_s

log_rank.R
 Calculating Log rank for all FDR genes

max_min.R
 Barchart and Overlap Barplot_2A.jpeg  Barplot_2B.jpeg  Barplot_S1.jpeg  Barplot_S2.jpeg

model_checks.R
 Evaluates the models generated from bootstrap method

model_checks_validation.R
 Validating  gene signature model results on CGGA data 

model_eval.R
 Plots for evaluating prognosis signature model

model.R
 Creates the SPG based gene signature model

overlap.R
 Calculating overlap with X, Y chromosome; gtex (tissue specific); sex-biased genes 

power_plot.R
 Correlation test and plots to estimate power of the experiment

prelim_data.R
 Exploratory data analysis of TCGA data

pval.R
 Histogram of p-values from cox-regression

Reports_all.ipynb
{

survival_ROC.R

