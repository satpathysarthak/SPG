#!/bin/bash
date
for i in {1..20}
do
   echo "$i"
   Rscript /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/20210425_model.R LGG $i 1>>/home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/log/20210512_model.log.out 2>> /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/log/20210512_model.log.err
   date
done
