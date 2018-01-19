#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 19/Dec/2016
# Last modified: 19/Dec/2015
# Optimised and customised to run at the Darwin HPC
# This is to make a symlink which is a dummy file name of different SLX

#SLX=SLX-11373 # 
#SLX_NEW=SLX-11367 # 
cd /home/ssg29/data/fastq/$SLX_NEW
for i in `ls ../$SLX/$SLX.NEB*.fq.gz`;do j=`echo $i | cut -d'/' -f3`; k="$SLX_NEW"${j/$SLX}; ln -s $i $k; done

#for i in `ls ../SLX-11371/SLX-11371.NEB*.fq.gz`;do j=`echo $i | cut -d'/' -f3`; k=SLX-11365${j/SLX-11371}; ln -s $i $k; done
#more ~/data/fastq/$SLX_NEW/$SLX.H*.contents.csv  | sed s/\"//g | awk 'BEGIN{FS=",";OFS=","}NR>1{print $2,$3,$4}' | sort
#more ~/data/fastq/$SLX_NEW/$SLX.H*.contents.csv  | sed s/\"//g | awk 'BEGIN{FS=",";OFS=","}NR>1{print $2,$3,$4}' | sort
