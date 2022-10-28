#!/bin/bash
date
for i in {61..80}
do
   echo "$i"
   Rscript /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/20210626_model.R HNSC $i 1>>/home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/log/20210626_model.log.out 2>> /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/log/20210626_model.log.err
   date
done
