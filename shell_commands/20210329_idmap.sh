#!/bin/sh
cd /home/workstation/Documents/Sarthak/gender_prognosis/cox_regression/FDR_gender_name/
for f in *_FDR_int.name
do
	date
	e=$( echo $f | sed 's/_FDR_int.name//')
	#script = $('python /home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap.py $e')
	#echo $script
	python /home/workstation/Documents/Sarthak/gender_prognosis/IDMapping/idmap.py $e
	date
	echo "$f and $e are presented"
done
