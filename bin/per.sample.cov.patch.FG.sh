#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 21/Oct/2015
# Last modified: 21/Oct/2015
# Optimised and customised to run at the Darwin HPC
# based on StringTie <https://github.com/gpertea/stringtie>

source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

PROJECT_DIR=$HOME/results/$SLX."Homo_sapiens.SE125.v2"

printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"

My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.bam # SE125 (JD)
#My_bam=$PROJECT_DIR/TopHat/$Barcode/final_merged_accepted_hits.bam # SE125 + SE50 (FG)
if [ ! -s $My_bam ];then 
	echo -e "\e[031m$My_bam not found. tophat failed? \e[0m\n"
	exit
fi

mkdir_unless $PROJECT_DIR/Coverage/$Barcode
# bam2coverage
#echo -e "bedtools genomecov $My_bam\n"
#time bedtools genomecov -split -ibam $My_bam > $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.txt
# bam2bedgraph
echo -e "bedtools genomecov -bg $My_bam\n"
time bedtools genomecov -bg -split -ibam $My_bam > $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.bedgraph
# make chrom.sizes
echo -e "samtools idxstats $My_bam\n"
time samtools idxstats $My_bam | awk 'BEGIN{OFS="\t"}!/*/{print $1,$2}' > $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.samtool.idxstats.chr.size.txt
# bedgraph2bw
echo -e "bedGraphToBigWig $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.bedgraph\n"
time bedGraphToBigWig $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.bedgraph $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.samtool.idxstats.chr.size.txt $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.bw

echo -e "\e[32m$Barcode done for StringTie & Bam2BigWig\e[0m" 
