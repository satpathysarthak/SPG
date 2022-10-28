#!/bin/sh
cd /home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_name/
for f in *_FDR_int.name
do
	date
	e=$( echo $f | sed 's/_FDR_int.name//')
	echo "Rscript /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/20210611_prelim_data.R $e"
	Rscript  /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/20210611_prelim_data.R $e
	date
	echo "$f and $e are presented"
done
