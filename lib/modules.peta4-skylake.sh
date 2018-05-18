#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 28/Feb/2018
# Last modified: 28/Feb/2018

# Python
#module load python-2.7.13-gcc-5.4.0-yubmrmn
module load python-3.6.1-gcc-5.4.0-xk7ym4l # for cutadapt --parallel

# Mapping
module load bowtie2-2.3.4.1-gcc-5.4.0-td2qanw
#module load samtools-1.4-gcc-5.4.0-derfxbk

# Methyl-Seq
#module load bismark/0.12.3

# RNA-Seq
#module add boost/1.55/boost_1.55.0-gcc-python_2.7.5
#module load tophat/2.0.12 # depends on boost above
#module load cufflinks/2.2.1

#bedtools 
#module load bedtools2-2.26.0-gcc-5.4.0-tqj36mo # ~/Install/bedtools

# Perl
#module load perl-5.24.1-gcc-5.4.0-dhq2pay

# R
#module load r-3.4.1-gcc-5.4.0-uj5r3tk # X11, jpeg, png disabled
module load r-3.4.1-gcc-5.4.0-jubrpyn # X11, jpeg, png enabled 


# Pandoc
module load pandoc/2.0.6
# TexLive for Rmarkdown
module load texlive/2015

# Java 
# for muTect - java 6 required
# for BisSNP - java 6 required
#module load java/jdk1.6.0_24

# for GATK - java>=7 required
# for fastqc
#module load java/jdk1.7.0_60
# java version "1.8.0_45" (peta4-skylake as of Feb/2018)

# ZLIB
module load zlib-1.2.11-gcc-5.4.0-dmjwhms

# Etc
module load htop
module load parallel 

