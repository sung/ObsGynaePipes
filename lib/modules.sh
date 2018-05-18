#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 5/Jul/2014
# Last modified: 8/Jul/2014 

# Python
module load python/2.7.5 # for cutadapt

# QC 
#module load fastqc/0.11.2		# moved to ~/Install/fastqc
#module load trim-galore/0.3.3  # moved to ~/Install/trim-galore
#module load cutadapt/1.4.2     # moved to ~/Install/cutadapt

# Mapping
module load bowtie/2.2.3
#module load samtools/0.1.19 # ~/Install/samtools
#module load gatk/3.1-1 # this will give $GATK_HOME (e.g. /usr/local/Cluster-Apps/gatk/3.1-1)
						# now installed at ~/Install/GenomeAnalysisTK-3.4-46-gbc02625/

# Methyl-Seq
module load bismark/0.12.3

# RNA-Seq
module add boost/1.55/boost_1.55.0-gcc-python_2.7.5
module load tophat/2.0.12 # depends on boost above
module load cufflinks/2.2.1

#bedtools 
module add bedtools/2.20.1

# Perl
module load perl/5.20.0 # default perl/5.10.1

# R
module load R/3.1.1 # before 3/Dec/2015 (Using Bioconductor version 3.0 (BiocInstaller 1.16.5), R version 3.1.1.)
					#  ggbio_1.14.0 ggplot2_2.1.0  data.table_1.10.4 gridExtra_0.9 (gridExtra_2.2.1)
#module load R/3.2.1 #  as of 3/Dec/2015 mainly for Bioconductor 3.1 (Bioconductor 3.2 works with R/3.2.2) library('cowplot')
					# as of 22/Jan/2016, upgraded to Bioconductor 3.2 via biocLite("BiocUpgrade")
					#  ggbio_1.18.5 ggplot2_2.2.1 data.table_1.10.4 gridExtra_2.2.1
#module load R/3.3.2 # 27/Apr/2017 (compatible with Bioconductor version 3.3 and 3.4)
					# Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.2 (2016-10-31).
#module load R/3.4.1-openblas # 11/Dec/2017 (Bioconductor version 3.6)
if [[ `hostname` =~ "login-e" ]]; then
    module load openblas-0.2.19-gcc-5.4.0-gstvz3w # this is to prevent libicuuc.so.42 from the peta4-skylake
    module load r-3.4.1-gcc-5.4.0-uj5r3tk # this is to prevent libicuuc.so.42 from the peta4-skylake
fi

# Pandoc
module load pandoc/1.17.2
# TexLive for Rmarkdown
#module load texlive/2015

# Etc
module load htop
module load parallel 

# Java 
# for muTect - java 6 required
# for BisSNP - java 6 required
#module load java/jdk1.6.0_24

# for GATK - java>=7 required
# for fastqc
module load java/jdk1.7.0_60
# java version "1.8.0_45" (peta4-skylake as of Feb/2018)


#htseq -- done ~/Install
#picard -- done ~/Install
#bissnp -- done ~/Install
#preseq
#bsExpress

# GCC G++
module load gcc/5.2.0 

# ZLIB
module load zlib/1.2.8
