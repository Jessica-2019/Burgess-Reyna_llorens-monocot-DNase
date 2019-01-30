#!/bin/bash

#This script runs pyDNAse wellington-bootstrap automatically 


#Set arguments

#################


echo -e "\n Enter the prefix of the treatment bam file\n"
read bam
echo -e "\n Enter the prefix of the control bam file\n"
read cbam


echo -e "\nPlease enter prefix of the treatment  1\n"
read t1
echo -e "\nPlease enter prefix of the control 1\n"
read c1



# Run PyDNASE for Tretment vs Control
echo -e "\n Running Wellington-bootstram Treatment vs Control\n"
wellington_bootstrap.py -fdrlimit -5 "$bam".bam "$cbam".bam "$t1".bed "$t1".TvsCdgf10 "$c1".TvsCdgf5


echo -e "\n Running Wellington-bootstram Control vs Treatment\n"
wellington_bootstrap.py -fdrlimit -5 "$cbam".bam "$bam".bam "$c1".bed  "$c1".CvsTdgf10 "$t1".CvsTdgf5 



echo -e "\n Running Wellington-bootstram Treatment vs Control\n"
wellington_bootstrap.py "$bam".bam "$cbam".bam "$t1".bed "$t1".TvsCdgf20 "$c1".TvsCdgf20 


echo -e "\n Running Wellington-bootstram Control vs Treatment\n"
wellington_bootstrap.py "$cbam".bam "$bam".bam "$c1".bed  "$c1".CvsTdgf20 "$t1".CvsTdgf20 



echo -e "\nwe are done\n"




wellington_bootstrap.py -fdrlimit -5 Bd_wl1.srt.bam Bd_nkd.srt.bam Bd_DHS.bed bd.TvsCdgf10 bd.TvsCdgf5


wellington_bootstrap.py -fdrlimit -5 "$bam".bam "$cbam".bam "$t1".bed "$t1".TvsCdgf10 "$c1".TvsCdgf5







