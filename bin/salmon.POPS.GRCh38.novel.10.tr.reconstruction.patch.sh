#!/bin/bash 

SPECIES='Homo_sapiens' # 'Homo_sapiens' or 'Mus_musculus' 
TR_PREFIX=GRCh38 # GRCh37|GRCh38
ENS_VER=82 # 82 (for manual downlaod)
SALMON_INDEX=$HOME/data/Salmon/transcriptome_index/$SPECIES/Ensembl/$TR_PREFIX/POPS.$TR_PREFIX.$ENS_VER.novel.10.tr.reconstruction.salmon.idx
SALMON_LIBTYPE=A # automatically infer the library type
GENOME=$HOME/data/genome/$SPECIES/Ensembl/$TR_PREFIX/Sequence/WholeGenomeFasta/genome.fa
GTF=/home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.gtf
GTF2=/home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.xloc.txt
GTF3=/home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.hgnc_symbol.txt
FASTA=/home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.fa
NT=12

#time gffread $GTF -g $GENOME -w $FASTA 
#time salmon index -t $GTF -i $SALMON_INDEX

SLX_IDS="SLX-9168 SLX-9169" # sand4
#SLX_IDS="SLX-10281 SLX-10402" # sand5
#SLX_IDS="SLX-9792 SLX-10284" # sand6
#SLX_IDS="SLX-10285 SLX-10283 SLX-10287" # sand7
#SLX_IDS="SLX-10402 SLX-10284 SLX-10287 SLX-9792 SLX-10281 SLX-10283 SLX-10285"
for SLX in $SLX_IDS; do
	PROJECT_DIR=$HOME/results/$SLX.$SPECIES.SE125.v2
	#PROJECT_DIR=$HOME/results/$SLX.$SPECIES.v1
	ALL_BARCODES=($(for i in `ls /scratch/obsgynae/POPS/data/fastq/$SLX/$SLX.*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f8 | cut -d'.' -f2 ; done | uniq))
	for Barcode in ${ALL_BARCODES[*]}; do
		echo -e "$SLX:$Barcode"
		SALMON_OUT=$PROJECT_DIR/Salmon/$Barcode; mkdir -p $SALMON_OUT
		TARGET=$SALMON_OUT/quant.genes.POPS.$TR_PREFIX.$ENS_VER.novel.10.tr.reconstruction.clean.txt
		if [ -s $TARGET ]; then
			awk '/TCONS/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.sf > $SALMON_OUT/$SLX.$Barcode.quant.POPS.$TR_PREFIX.$ENS_VER.novel.10.tr.reconstruction.clean.txt
			mv $TARGET $SALMON_OUT/$SLX.$Barcode.quant.genes.POPS.$TR_PREFIX.$ENS_VER.novel.10.tr.reconstruction.clean.txt
		else
			FastQ_array=() #initialise
			for FastQ_file in `ls /scratch/obsgynae/POPS/data/fastq/$SLX/$SLX.$Barcode*.fq.gz | grep -v lost`; do FastQ_array+=($FastQ_file); done
			Merged_FastQ=$(printf " %s" "${FastQ_array[@]}")
			Merged_FastQ=${Merged_FastQ:1} # remove the first space 
			echo -e "salmon quant -p $NT -i $SALMON_INDEX -l $SALMON_LIBTYPE -r $Merged_FastQ -g $GTF2 -o $SALMON_OUT"
			time salmon quant -p $NT -i $SALMON_INDEX -l $SALMON_LIBTYPE -r $Merged_FastQ -g $GTF2 -o $SALMON_OUT
			awk '/TCONS/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.sf > $SALMON_OUT/$SLX.$Barcode.quant.POPS.$TR_PREFIX.$ENS_VER.novel.10.tr.reconstruction.clean.txt
			awk '/XLOC/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.genes.sf > $SALMON_OUT/$SLX.$Barcode.quant.genes.POPS.$TR_PREFIX.$ENS_VER.novel.10.tr.reconstruction.clean.txt
		fi
	done
done
