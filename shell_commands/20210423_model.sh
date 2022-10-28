#!/bin/bash
for i in {1..1000}
do
   echo "$i"
   Rscript /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/20210423_model.R GBM $i 1>>/home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/log/20210423_model.log.out 2>> /home/workstation/Documents/Sarthak/gender_prognosis/Rcodes/log/20210423_model.log.err
done
