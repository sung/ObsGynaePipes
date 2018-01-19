#!/bin/bash
#CSMD1
MY_GENE=CSMD1
REF=GRCh37 #GRCh37|GRCh38
ENS=75     #75|82
OUT_DIR=/home/ssg29/data/genome/Homo_sapiens/Ensembl/$REF/Annotation/Genes/
GTF=$OUT_DIR/Homo_sapiens.$REF.$ENS.gtf
UNIQ_GTF=$OUT_DIR/Homo_sapiens.$REF.$ENS.$MY_GENE.exon.union.sorted.gtf
MERGED_BED=$OUT_DIR/Homo_sapiens.$REF.$ENS.$MY_GENE.exon.union.sorted.reduced.bed

# unique exon region (uniq start and end position)
time grep $MY_GENE $GTF | awk '$3=="exon"{print $_}'| \
	sort -k1,1 -k4,4n -k5,5n | \
	# group by 'chr', 'start', 'end'
	bedtools groupby -grp 1,4,5 -opCols 6 -ops first -full | \
	awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > $UNIQ_GTF

##########################################
## Merge (or Reduce) unique exon region ##
##########################################
bedtools merge -s -c 9 -o collapse -delim "; " -i $UNIQ_GTF | \
	awk 'BEGIN{FS="\t";OFS="\t"}{split($4,a,";"); ense=""; for(i in a) if(a[i]~/exon_id/){split(a[i], b, "exon_id"); ense=ense""b[2]} print $1,$2,$3,"exon_id="ense}' | \
	sed 's/"//g' > $MERGED_GTF
