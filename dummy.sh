#!/bin/bash

FLAG=$(stat -c%s "/home/ssg29/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-9169.Homo_sapiens.SE125.v2/BWA/D712_D508/SLX-9169.D712_D508.C80CGANXX.s_8.bwa.sam") # file size
if [ $FLAG -lt 1200 ];then
    echo -e "skip" 
else
    echo -e "run" 
fi

FILE="/home/ssg29/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10283.Homo_sapiens.v1/CIRI2/D710_D507/SLX-10283.D710_D507.merged.bwa.sam.CIRI2.txt"
if [ -s $FILE ]; then
    echo -e "yes"
else
    echo -e "no"
fi
