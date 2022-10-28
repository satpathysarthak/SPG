#!/bin/sh
cd /home/workstation/Documents/Sarthak/gender_prognosis/analysis/bivariate_go/reports 
for f in *_report.xlsx
do
	date
	e=$( echo $f | sed 's/_report.xlsx//')
	#script = $('python /home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap.py $e')
	#echo $script
	python /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/20210406_go_report.py $e
	date
	echo "$f and $e are presented"
done
