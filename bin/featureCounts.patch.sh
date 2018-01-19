#!/bin/bash
# Author: Sung Gong <sung@bio.cc>
# Last modified in 25/Oct/2016
source $HOME/lib/sung.sh # mkdir_unless

REF="GRCh37" # GRCh38
ENS=75 #82
ENS_VER=$REF.$ENS
NT=4
MY_GENE="CSMD1"
GTF=$HOME/data/genome/Homo_sapiens/Ensembl/$REF/Annotation/Genes/Homo_sapiens.$ENS_VER.$MY_GENE.exon.union.sorted.gtf # see bin/fetch.exon.union.sh
#GTF=$HOME/results/RNA-Seq/Placentome/Cuffcompare/CTL/CTL.$ENS_VER.cuffcompare.$MY_GENE.XLOC_069068.exon.union.sorted.gtf # based on the reconstructed Placenta transcriptome
SOURCE="JD" #FG

#JD
#ALL_SLX="SLX-10281 SLX-10402 SLX-9792 SLX-10284 SLX-10285 SLX-10283 SLX-10287"
#ALL_SLX="SLX-10281 SLX-10402 SLX-9792 SLX-10284"
#ALL_SLX="SLX-10285 SLX-10283 SLX-10287"

#FG
#ALL_SLX="SLX-9168 SLX-9169"
#ALL_SLX="SLX-9168"
for SLX in $ALL_SLX
do
	if [[ $SOURCE == "FG" ]]; then
		if [[ $REF == "GRCh37" ]]; then
			SUBREAD_DIR=$HOME/results/$SLX.Homo_sapiens.SE125.v1/subRead #FG
			TOPHAT_DIR=$HOME/results/$SLX.Homo_sapiens.SE125.v1/TopHat   #FG
		else
			SUBREAD_DIR=$HOME/results/$SLX.Homo_sapiens.SE125.v2/subRead #FG
			TOPHAT_DIR=$HOME/results/$SLX.Homo_sapiens.SE125.v2/TopHat   #FG
		fi
	else 
		if [[ $REF == "GRCh37" ]]; then
			SUBREAD_DIR=$HOME/results/$SLX.Homo_sapiens.v2/subRead        #JD
			TOPHAT_DIR=$HOME/results/$SLX.Homo_sapiens.v2/TopHat          #JD
		else
			SUBREAD_DIR=$HOME/results/$SLX.Homo_sapiens.v1/subRead        #JD
			TOPHAT_DIR=$HOME/results/$SLX.Homo_sapiens.v1/TopHat          #JD
		fi
	fi

	mkdir_unless $SUBREAD_DIR
	BARCODE=$(for i in `ls /home/ssg29/scratch/data/fastq/$SLX/$SLX*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f8 | cut -d'.' -f2 ; done | uniq)
	for i in $BARCODE
	do 
		echo $SLX"."$i
		mkdir_unless $SUBREAD_DIR/$i
		MY_CNT_FILE=$SUBREAD_DIR/$i/$SLX.$i.featureCounts.$ENS_VER.$MY_GENE.exon_id.txt
		#MY_CNT_FILE=$SUBREAD_DIR/$i/$SLX.$i.featureCounts.$ENS_VER.$MY_GENE.XLOC_069068.exon_id.txt # Placentome (GRCh38) 

		# -T: No. of core
		# -a: GTF file
		# -s: 2 (reversely stranded)
		# -f: Perform read counting at feature level (eg. counting reads for exons rather than genes).
		# -O: Assign reads to all their overlapping meta-features (or features if -f is specified).
		echo -e "featureCounts -T $NT -g exon_id -a $GTF -s 2 -f -O\n"
		time featureCounts -T $NT \
			-g exon_id \
			-a $GTF \
			-s 2 \
			-f \
			-O \
			-o $MY_CNT_FILE \
			$TOPHAT_DIR/$i/merged_accepted_hits.bam
		time awk 'BEGIN{OFS="\t"}NR>2{print $1,$7}' $MY_CNT_FILE | sort > ${MY_CNT_FILE%.txt}.clean.txt
	done
done
