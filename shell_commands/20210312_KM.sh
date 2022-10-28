#!/bin/bash
while read -r line; do
	# Reading each line
	echo $line
	Rscript /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/20201025_KM_plots-copy.R $line 
done < /home/workstation/Documents/Sarthak/gender_prognosis/analysis/trial.txt
#done < /home/workstation/Documents/Sarthak/gender_prognosis/analysis/max_min_reg.txt
