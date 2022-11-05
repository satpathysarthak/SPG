# SPG
Scripts and Processed Data for SPG manuscript
==> 20201027_list_gender.R <==
## Script to generate plots for Gene Expression

==> 20201028_KM_genes.R <==
## KM plot superstar genes expession

==> 20201112_pval.R <==
input_dir <- "/media/sarthak/gender_prognosis/cox_regression/gender"

==> 20201124_FDR_gender.R <==
## FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.20

==> 20210131_power_plot.R <==
setwd("/home/workstation/Documents/Sarthak/gender_prognosis")

==> 20210213_cox_gender_z.R <==
## Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated

==> 20210215_FDR_gender_zscore.R <==
## FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.20

==> 20210215_gender_pub_plots.R <==
# Barchart and Overlap Barplot_2A.jpeg  Barplot_2B.jpeg  Barplot_S1.jpeg  Barplot_S2.jpeg

==> 20210220_go_plots.R <==
library(topGO)

==> 20210220_max_min.R <==
# Barchart and Overlap Barplot_2A.jpeg  Barplot_2B.jpeg  Barplot_S1.jpeg  Barplot_S2.jpeg

==> 20210221_idmap_KM.R <==
#ID mapping: Specific from taking input from .int file

==> 20210309_int_est_overlap.R <==
setwd("/home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_zscore")

==> 20210310_heatmap.R <==
## Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated

==> 20210312_log_rank.R <==
#### Calculating Log rank for all FDR genes

==> 20210313_log_rank_table.R <==
library(data.table)

==> 20210317_FDR_bivariate.R <==
## Script to perform Survival Analysis with Tumor data using 2 models: Male Female

==> 20210318_FDR_bivariate_table.R <==
library(data.table)

==> 20210325_FDR_bivariate_plots.R <==
library(data.table)

==> 20210331_FDR_bivariate_GO_create.R <==
output_dir <- "/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/"

==> 20210406_go_report.py <==
import pandas as pd

==> 20210408_bivariate_GO_plots.R <==
dir_path <- '/home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/reports_all'

==> 20210413_int_beta_scatter.R <==
library(data.table)

==> 20210423_model.R <==
rm(list = ls())

==> 20210424_model.R <==
rm(list = ls())

==> 20210425_model.R <==
rm(list = ls())

==> 20210512_model_checks.r <==
rm(list = ls())

==> 20210518_model_checks_validation.r <==
### Validating ressssssults on CGGA data 

==> 20210529_HNSC_data_processing.R <==
## Microarray data processing

==> 20210530_HNSC_validation_KM.R <==
# plotting KM  plots

==> 20210611_prelim_data.R <==
##Barchart for number of significant genes with p<0.20 FDR cutoff

==> 20210618_model_eval.R <==
## baplots

==> 20210623_model.R <==
rm(list = ls())

==> 20210626_model.R <==
rm(list = ls())

==> 20220523_heatmap_IDH.R <==
## Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated

==> cox_gender_miR.R <==
##Script automated for miR gender Cox proportional regression

==> cox_gender.R <==
## Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated

==> cox_menopausal.R <==
## Script to perform Survival Analysis with Tumor data using 2 models: Interaction and Estimated

==> cox_uni_gender.R <==
## Script to perform Survival Analysis with Tumor data using univariate models on subdata

==> create_processed_sample.R <==
## Merge information on gender with Clinical-OS: Creating processed Clinical Sample File

==> FDR_gender_miR.R <==
## miRNA regression: FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 

==> FDR_plots.R <==
## Plot for miRNA and mRNA expreesion: FDR Filter to get regression results which qualify the mean absolute cutoff and FDR

==> FDR_plots_zscore.R <==
## FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.10. Plot the changes. 

==> FDR_uni_gender.R <==
## FDR Filter to get regression results which qualify the mean absolute cutoff and FDR 0.20

==> gender_pub_plots.R <==
# Barchart and Overlap Barplot_2A.jpeg  Barplot_2B.jpeg  Barplot_S1.jpeg  Barplot_S2.jpeg

==> go_rnk_file_creation.R <==
#ID mapping: Specific from taking input from .int file

==> gsea_rnk_file_creation.R <==
#ID mapping: Specific from taking input from .int file

==> hist_sample.R <==
## Script to generate plots for Gene Expression

==> HNSC_validation_AUC.R <==
### Calculates AUC for HNSC validation data

==> idmap_fun.RData <==
�      �X�k�F?����6Q�lӼ���Ą�2J�H�<B�P���0�-'Jl�H2Y	���c��Ż�=��'y�o��x_�����{z��O��'�'�Iʒ䯧͝�	I(T���2������U�)��8F�B�P)���05��`h��4��3�e���-�M�p)��*���͎�Yϴ���-�r��y�z��:V��n{8�������7��V�-�z��F�qZ���ʵ�*�0���k1[􌭞cR�~i����\�^�Y�ڝ���<�~�v���u�°���:��e��[?:��z`Zu������>/�my��Y�ɰ��1=n`�`����v��9�D���_�y�i4~>��F�Et�f��َ�7�f5x�#ˤr�����l�rEϣｶ��B�~D�!{L���v���h衠�B����0uN]���/*�ų�9}��R'a�5�˳g �ɉC�����PL��{Q��P����ȱ�k��U�5v~�`����;>foY��v�[��2��c?�/,@�ǂ�IE�cOpة�7_���:�!��:����H��}q�/1^ī	���>N��dE5�Ƶ@�Ҩ	@H�	>$Q����;�߁��FM�� ��QӀ�Dx�i��D]F��?yE�) 1Z��Si�U@z.<�U��D]F_�

==> idmap_one.R <==
idmap <- function(name_ensembl){

==> idmap_table.R <==
input_file_ref_from <- "/home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap_ensembl.tab"

==> KM_plots.R <==
# Script for plotting KM Plots`

==> overlap.R <==
## Calculating overlap with X, Y chromosome; gtex (tissue specific); sex-biased genes 

==> Reports_all.ipynb <==
{

==> survival_ROC.R <==

