#!/bin/bash

#This script runs pyDNAse wellington automatically 


#Set arguments

#################



echo -e "\n Enter the prefix of the bam1 file\n"
read bam1
echo -e "\n Enter the prefix of the bam2 file\n"
read bam2
echo -e "\n Enter the prefix of the bam3 file\n"
read bam3
echo -e "\n Enter the prefix of the bam4 file\n"
read bam4
echo -e "\n Enter the prefix of the bam5 file\n"
read bam5
echo -e "\n Enter the prefix of the bam6 file\n"
read bam6



echo -e "\nPlease enter prefix of the bed file 1  1\n"
read t1
echo -e "\nPlease enter prefix of the bed file 2  1\n"
read t2
echo -e "\nPlease enter prefix of the bed file 3  1\n"
read t3
echo -e "\nPlease enter prefix of the bed file 4  1\n"
read t4
echo -e "\nPlease enter prefix of the bed file 5  1\n"
read t5
echo -e "\nPlease enter prefix of the bed file 6  1\n"
read t6



# Run wellington
 
mkdir ./"$t1".dgf
echo -e "\n Running Wellington-bootstram 1 \n"
wellington_footprints.py -fdrlimit -5 "$t1".bed "$bam1".bam "$t1".dgf

mkdir ./"$t2".dgf
echo -e "\n Running Wellington-bootstram 2 \n"
wellington_footprints.py -fdrlimit -5 "$t2".bed "$bam2".bam "$t2".dgf 

mkdir ./"$t3".dgf
echo -e "\n Running Wellington-bootstram 3 \n"
wellington_footprints.py -fdrlimit -5  "$t3".bed "$bam3".bam "$t3".dgf 


mkdir ./"$t4".dgf
echo -e "\n Running Wellington-bootstram 4 \n"
wellington_footprints.py -fdrlimit -5 "$t4".bed "$bam4".bam "$t4".dgf 

mkdir ./"$t5".dgf
echo -e "\n Running Wellington-bootstram 5 \n"
wellington_footprints.py -fdrlimit -5 "$t5".bed "$bam5".bam "$t5".dgf 

mkdir ./"$t6".dgf
echo -e "\n Running Wellington-bootstram 6 \n"
wellington_footprints.py -fdrlimit -5 "$t6".bed "$bam6".bam "$t6".dgf 

echo -e "\nwe are done\n"














