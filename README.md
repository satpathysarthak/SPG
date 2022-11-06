# SPG

The repository contains Scripts and Processed Data for SPG manuscript.

## Cox Regression model results

cox_regression/ has all the cox regression models.
cox_regression/FDR_gender_zscore/ has the cox-regression model which 

## Scripts


bivariate_GO_plots.R: Plots from processed GO analysis

cox_gender_miR.R: Script automated for miR gender Cox proportional regression

cox_gender.R: Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated

cox_gender_z.R: Script to perform Survival Analysis with Tumor data post z-score transformation using 2 models: Interaction and Estimated

cox_menopausal.R: Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated

cox_uni_gender.R: Script to perform Survival Analysis with Tumor data using univariate models on subdata

create_processed_sample.R: Merge information on gender with Clinical-OS: Creating processed Clinical Sample File

FDR_bivariate_GO_create.R: Creates a ranked file from output of FunrichR

FDR_bivariate_plots.R: Plot female vs male beta for ns_ns,s_ns,ns_s,s_s

FDR_bivariate.R: Script to perform Survival Analysis with Tumor data using 2 models: Male Female

FDR_bivariate_table.R: Genes in the category: ns_ns,s_ns,ns_s,s_s

FDR_gender_miR.R: miRNA regression: FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 

FDR_gender.R:  FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.10

FDR_gender_zscore.R: FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.10

FDR_plots.R: Plot for miRNA and mRNA expreesion: FDR Filter to get regression results which qualify the mean absolute cutoff and FDR

FDR_plots_zscore.R: FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.10. Plot the changes. 

FDR_uni_gender.R: FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.20

gender_pub_plots.R: Barchart and Overlap Barplot_2A.jpeg  Barplot_2B.jpeg  Barplot_S1.jpeg  Barplot_S2.jpeg

go_plots.R: Plots from Enrichment scores

go_report.py: Processing the reports generated by FunrichR

heatmap_IDH.R: Heatmap for gene expression betwwen IDH status

heatmap.R:nHeatmap of gene expression for genes from Survial Analysis

hist_sample.R: Script to generate plots for Gene Expression

HNSC_data_processing.R: Microarray data processing

HNSC_validation_AUC.R:  Calculates AUC for HNSC validation data

HNSC_validation_KM.R:  KM  plots for HNSC validation and training data


idmap_KM.R: ID mapping: Specific from taking input from .int file

int_beta_scatter.R: Creates a scatter plot male and female coefficients 

int_est_overlap.R: Calculating overlap between Main-effect and Interaction model 

KM_genes.R:  KM plot specific genes

KM_plots.R: Script for plotting KM Plots`

list_gender.R: Script to generate plots for data exploration: ditribution for cencored and gender data`

ListOfSigNS.R: Counting ns_ns,s_ns,ns_s,s_s

log_rank.R: Calculating Log rank for all FDR genes

max_min.R: Barchart and Overlap Barplot_2A.jpeg  Barplot_2B.jpeg  Barplot_S1.jpeg  Barplot_S2.jpeg

model_checks.R: Evaluates the models generated from bootstrap method

model_checks_validation.R: Validating  gene signature model results on CGGA data 

model_eval.R: Plots for evaluating prognosis signature model

model.R: Creates the SPG based gene signature model

overlap.R: Calculating overlap with X, Y chromosome; gtex (tissue specific); sex-biased genes 

power_plot.R: Correlation test and plots to estimate power of the experiment

prelim_data.R: Exploratory data analysis of TCGA data

pval.R" Histogram of p-values from cox-regression

