#!/bin/bash
date
for i in {1..100}
do
   echo "$i"
   Rscript /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/20210424_model.R GBM $i 1>>/home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/log/20210424_model.log.out 2>> /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/log/20210424_model.log.err
   date
done
