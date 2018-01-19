#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 23/Oct/2015
# Last modified: 23/Oct/2015
# Optimised and customised to run at the Darwin HPC

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

PROJECT=Placentome
PROJECT_DIR=$HOME/results/RNA-Seq/$PROJECT

COHORT="MY_COHORT" # CTLF, CTLM, POPS, CTL, SGA, PET
echo -e "COHORT=$COHORT"

mkdir_unless $PROJECT_DIR/Cuffcompare
mkdir_unless $PROJECT_DIR/Cuffcompare/$COHORT

################################################### 
## READ ~/results/RNA-Seq/Placentome/Meta/README ##
################################################### 

########
## FG ##
########
#awk 'BEGIN{FS=","}$3==0{print "/home/ssg29/results/"$1".Homo_sapiens.SE125.v1/StringTie/"$2"/"$1".Homo_sapiens.SE125.v1."$2".stringtie.GRCh37.75.assembly.gtf"}' $PROJECT_DIR/Meta/FG.barcode.txt > $PROJECT_DIR/Meta/FG.SE125.SE50.$COHORT.stringtie.GRCh37.75.gtf.list 
#awk 'BEGIN{FS=","}$3==0{print "/home/ssg29/results/"$1".Homo_sapiens.SE125.v2/StringTie/"$2"/"$1".Homo_sapiens.SE125.v2."$2".stringtie.GRCh38.82.assembly.gtf"}' $PROJECT_DIR/Meta/FG.barcode.txt > $PROJECT_DIR/Meta/FG.SE125.SE50.$COHORT.stringtie.GRCh38.82.gtf.list 

########
## JD ##
########
#awk 'BEGIN{FS=","}$3!~/P/{print "/home/ssg29/results/"$1".Homo_sapiens.v1/StringTie/"$2"/"$1".Homo_sapiens.v1."$2".stringtie.GRCh38.82.assembly.gtf"}' $PROJECT_DIR/Meta/JD.barcode.txt > $PROJECT_DIR/Meta/JD.$COHORT.stringtie.GRCh38.82.gtf.list 

################
## Placentome ##
################
#cat ~/results/RNA-Seq/SGA.AGA/Meta/FG.SE125.SE50.AGA.stringtie.GRCh38.82.gtf.list ~/results/RNA-Seq/PET/Meta/JD.AGA.stringtie.GRCh38.82.gtf.list > ~/results/RNA-Seq/Placentome/Meta/CTL.stringtie.GRCh38.82.gtf.list

####################
## Healthy Female ##
####################
#awk 'BEGIN{FS=","}NR>1 && $6==0 && $1!="84C" && $1!="88" && $1!="80C" && $1!="93C" && $1!="94C" && $7=="F"{split($2,a,"/"); myFile="/home/ssg29/results/"a[5]"/StringTie/"a[7]"/"a[5]"."a[7]".stringtie.GRCh37.75.assembly.gtf"; print myFile}' ~/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/Meta/meta.ALL.csv > ~/results/RNA-Seq/Placentome/Meta/CTLF.stringtie.GRCh37.75.gtf.list

##################
## Healthy Male ##
##################
#awk 'BEGIN{FS=","}NR>1 && $6==0 && $1!="84C" && $1!="88" && $1!="80C" && $1!="93C" && $1!="94C" && $7=="M"{split($2,a,"/"); myFile="/home/ssg29/results/"a[5]"/StringTie/"a[7]"/"a[5]"."a[7]".stringtie.GRCh37.75.assembly.gtf"; print myFile}' ~/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/Meta/meta.ALL.csv > ~/results/RNA-Seq/Placentome/Meta/CTLM.stringtie.GRCh37.75.gtf.list

if [ -s $PROJECT_DIR/Meta/$COHORT.stringtie.$TR_PREFIX.$ENS_VER.gtf.list ];then
	# cuffcompare v2.2.1 (4237)
	# -C: include the "contained" transcripts in the .combined.gtf file
	# -G: generic GFF input file(s): do not assume Cufflinks GTF, do not
	#     discard any intron-redundant transfrags)
	# -T: do not generate .tmap and .refmap files for each input file
	# -V: verbose processing mode (showing all GFF parsing warnings)
	echo -e "cuffcompare -i $PROJECT_DIR/Meta/$COHORT.stringtie.$TR_PREFIX.$ENS_VER.gtf.list\n"
	time cuffcompare \
		-i $PROJECT_DIR/Meta/$COHORT.stringtie.$TR_PREFIX.$ENS_VER.gtf.list \
		-r $GTF \
		-s $GENOME \
		-C \
		-G \
		-T \
		-V \
		-o $PROJECT_DIR/Cuffcompare/$COHORT/$COHORT.$TR_PREFIX.$ENS_VER.cuffcompare \
		2> $PROJECT_DIR/Cuffcompare/$COHORT/$COHORT.$TR_PREFIX.$ENS_VER.cuffcompare.log
else
	echo -e "\e[031$PROJECT_DIR/Meta/$COHORT.stringtie.$TR_PREFIX.$ENS_VER.gtf.list not found\e[0m\n"
	exit
fi
echo -e "Cuffcompare done"

My_GTF=$PROJECT_DIR/Cuffcompare/$COHORT/$COHORT.$TR_PREFIX.$ENS_VER.cuffcompare.combined.gtf
if [ -s $My_GTF ];then
	echo -e "indexing $My_GTF"
	(grep ^"#" $My_GTF; grep -v ^"#" $My_GTF | sort -k1,1 -k4,4n) | bgzip > ${My_GTF%.gtf}.sorted.gtf.gz

	# https://github.com/dasmoth/gtf2bed
	echo -e "GTF2BED"
	echo -e "java -jar $HOME/Install/gtf2bed/target/gtf2bed-0.0.4-SNAPSHOT-standalone.jar $My_GTF > ${My_GTF%.gtf}.geneName.bed\n"
	time java -jar $HOME/Install/gtf2bed/target/gtf2bed-0.0.4-SNAPSHOT-standalone.jar $My_GTF > ${My_GTF%.gtf}.geneName.bed

fi

if [ -s ${My_GTF%.gtf}.geneName.bed ];then
	# replace GeneName(v14) with GeneId(v13)
	time awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$13,$15,$16,$17}' ${My_GTF%.gtf}.geneName.bed > ${My_GTF%.gtf}.bed

	#http://hgdownload.cse.ucsc.edu/admin/exe/
	echo -e "Bed2BigBed"
	time bedToBigBed -type=bed12+4 -as=$HOME/Install/gtf2bed/gencode.as -extraIndex=name ${My_GTF%.gtf}.geneName.bed $CHRINFO ${My_GTF%.gtf}.geneName.bb 
	time bedToBigBed -type=bed12+4 -as=$HOME/Install/gtf2bed/gencode.as -extraIndex=name ${My_GTF%.gtf}.bed $CHRINFO ${My_GTF%.gtf}.bb 

	echo -e "Gene Name"
	time cut -f 4,14 ${My_GTF%.gtf}.geneName.bed > ${My_GTF%.gtf}.geneName.txt
	# index it
	time ixIxx ${My_GTF%.gtf}.geneName.txt ${My_GTF%.gtf}.geneName.ix ${My_GTF%.gtf}.geneName.ixx 
fi

echo -e "Cuffcompare done for $COHORT"
