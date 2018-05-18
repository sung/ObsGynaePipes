#!/bin/bash
# Author: Sung Gong <sung@bio.cc>
# Last modified in 25/Oct/2016

REF="GRCh38" # GRCh38, GRCh37
ENS=82 #82, 75
ENS_VER=$REF.$ENS
NT=16
MY_GENE="CSMD1"
#GTF=$HOME/data/genome/Homo_sapiens/Ensembl/$REF/Annotation/Genes/Homo_sapiens.$ENS_VER.$MY_GENE.exon.union.sorted.gtf # see bin/fetch.exon.union.sh
#GTF=$HOME/results/RNA-Seq/Placentome/Cuffcompare/CTL/CTL.$ENS_VER.cuffcompare.$MY_GENE.XLOC_069068.exon.union.sorted.gtf # based on the reconstructed Placenta transcriptome
GTF=$HOME/data/genome/Homo_sapiens/Ensembl/$REF/Annotation/Genes/Homo_sapiens.$ENS_VER.gtf # Homo_sapiens.GRCh37.82.UCSC.gtf 
SOURCE="JD" #FG

#JD
#ALL_SLX="SLX-10281 SLX-10402 SLX-9792 SLX-10284 SLX-10285 SLX-10283 SLX-10287"
ALL_SLX="SLX-10281 SLX-10402 SLX-9792 SLX-10284"
#ALL_SLX="SLX-10285 SLX-10283 SLX-10287"

#FG
#ALL_SLX="SLX-9168 SLX-9169"
#ALL_SLX="SLX-9168"
for SLX in $ALL_SLX
do
	if [[ $SOURCE == "FG" ]]; then
		if [[ $REF == "GRCh37" ]]; then
			FTCNT_DIR=$HOME/results/$SLX.Homo_sapiens.SE125.v1/featureCount #FG
			TOPHAT_DIR=$HOME/results/$SLX.Homo_sapiens.SE125.v1/TopHat   #FG
		else # GRCh38
			FTCNT_DIR=$HOME/results/$SLX.Homo_sapiens.SE125.v2/featureCount #FG
			TOPHAT_DIR=$HOME/results/$SLX.Homo_sapiens.SE125.v2/TopHat   #FG
		fi
	else 
		if [[ $REF == "GRCh37" ]]; then
			FTCNT_DIR=$HOME/results/$SLX.Homo_sapiens.v2/featureCount        #JD
			TOPHAT_DIR=$HOME/results/$SLX.Homo_sapiens.v2/TopHat          #JD
		else # GRCh38
			FTCNT_DIR=$HOME/results/$SLX.Homo_sapiens.v1/featureCount        #JD
			TOPHAT_DIR=$HOME/results/$SLX.Homo_sapiens.v1/TopHat          #JD
		fi
	fi

	mkdir -p $FTCNT_DIR
	BARCODE=$(for i in `ls /home/ssg29/data/fastq/$SLX/$SLX*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f7 | cut -d'.' -f2 ; done | uniq)
	for i in $BARCODE
	do 
		echo $SLX"."$i
		mkdir -p $FTCNT_DIR/$i
		#MY_CNT_FILE=$FTCNT_DIR/$i/$SLX.$i.featureCounts.$ENS_VER.$MY_GENE.exon_id.txt
		#MY_CNT_FILE=$FTCNT_DIR/$i/$SLX.$i.featureCounts.$ENS_VER.$MY_GENE.XLOC_069068.exon_id.txt # Placentome (GRCh38) 
        MY_CNT_FILE=$FTCNT_DIR/$i/$SLX.$i.featureCount.unmap.rescued.$ENS_VER.txt

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
