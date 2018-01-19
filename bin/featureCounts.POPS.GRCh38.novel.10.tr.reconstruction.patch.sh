#!/bin/bash 

#GTF="/home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.gtf"
GTF="/home/ssg29/data/genome/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.82.gtf"
NT=8

# featureCounts
# -s
# 1: for SMARTer (clone-tech plasma)
# 2: for placenta tissue assay (TrueSeq total RNA-Seq)

#featureCounts -T $NT -a $GTF -Q 10 -s 2 -p -C -o  ~/results/SLX-9168.Homo_sapiens.SE125.v2/featureCount/D701_D502/SLX-9168.D701_D502.POPS.GRCh38.novel.10.tr.reconstruction.txt  ~/results/SLX-9168.Homo_sapiens.SE125.v2/TopHat/D701_D502/merged_accepted_hits.bam
#time awk 'BEGIN{OFS="\t"}NR>2{print $1,$7}' ~/results/SLX-9168.Homo_sapiens.SE125.v2/featureCount/D701_D502/SLX-9168.D701_D502.POPS.GRCh38.novel.10.tr.reconstruction.txt | sort > ~/results/SLX-9168.Homo_sapiens.SE125.v2/featureCount/D701_D502/SLX-9168.D701_D502.POPS.GRCh38.novel.10.tr.reconstruction.clean.txt

# GTF to FASTA
#time gffread /home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.gtf -g ~/data/genome/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa -w /home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.fa

#SLX_IDS="SLX-9168 SLX-9169"
#for SLX in $SLX_IDS;do
#	echo -e $SLX
#	#for i in `ls ~/results/$SLX.Homo_sapiens.SE125.v2/TopHat/*/merged_accepted_hits.bam`; do echo $i; j=`echo $i | cut -d'/' -f7`; k=~/results/$SLX.Homo_sapiens.SE125.v2/featureCount/$j; FILE=$k/$SLX.$j.POPS.GRCh38.novel.10.tr.reconstruction.txt; mkdir -p $k; time featureCounts -T $NT -a $GTF -Q 10 -s 2 -p -C -o $FILE $i; time awk 'BEGIN{OFS="\t"}NR>2{print $1,$7}' $FILE | sort > ${FILE%txt}clean.txt ; done
#	for i in `ls ~/results/$SLX.Homo_sapiens.SE125.v2/TopHat/*/merged_accepted_hits.bam`; do echo $i; j=`echo $i | cut -d'/' -f7`; k=~/results/$SLX.Homo_sapiens.SE125.v2/featureCount/$j; FILE=$k/$SLX.$j.featureCount.unmapped.rescue.txt; mkdir -p $k; time featureCounts -T $NT -a $GTF -Q 10 -s 2 -o $FILE $i; time awk 'BEGIN{OFS="\t"}NR>2{print $1,$7}' $FILE | sort -k1,1 > ${FILE%.txt}.clean.txt ; done
#done

SLX_IDS="SLX-10281 SLX-10402 SLX-9792 SLX-10284 SLX-10285 SLX-10283 SLX-10287"
for SLX in $SLX_IDS;do
	echo -e $SLX
	#for i in `ls ~/results/$SLX.Homo_sapiens.v1/TopHat/*/merged_accepted_hits.bam`; do echo $i; j=`echo $i | cut -d'/' -f7`; k=~/results/$SLX.Homo_sapiens.v1/featureCount/$j; FILE=$k/$SLX.$j.POPS.GRCh38.novel.10.tr.reconstruction.txt; mkdir -p $k; time featureCounts -T $NT -a $GTF -Q 10 -s 2 -p -C -o $FILE $i; time awk 'BEGIN{OFS="\t"}NR>2{print $1,$7}' $FILE | sort > ${FILE%txt}clean.txt ; done
	for i in `ls ~/results/$SLX.Homo_sapiens.v1/TopHat/*/merged_accepted_hits.bam`; do echo $i; j=`echo $i | cut -d'/' -f7`; k=~/results/$SLX.Homo_sapiens.v1/featureCount/$j; FILE=$k/$SLX.$j.featureCount.unmapped.rescue.txt; mkdir -p $k; time featureCounts -T $NT -a $GTF -Q 10 -s 2 -o $FILE $i; time awk 'BEGIN{OFS="\t"}NR>2{print $1,$7}' $FILE | sort -k1,1 > ${FILE%.txt}.clean.txt ; done
done
