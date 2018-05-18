#!/bin/bash

    NT=16 
    SPECIES='Homo_sapiens' # 'Homo_sapiens' or 'Mus_musculus' or 'PhiX'
    SLX="SLX-Xten2017"
    Cell="XXX"
    Barcode="23409_1"
    BWA_INDEX_BASE=$HOME/data/genome/$SPECIES/Ensembl/GRCh38/Sequence/BWAIndex/genome # initial mapping

	Merged_Fw_FastQ=$HOME/results/$SLX.$SPECIES.v1/Trim/$Barcode/$SLX.$Barcode.1K.s_x.r_1.fq.gz
	Merged_Rv_FastQ=$HOME/results/$SLX.$SPECIES.v1/Trim/$Barcode/$SLX.$Barcode.1K.s_x.r_2.fq.gz
	TRIMMED_FASTQ_FILES="$Merged_Fw_FastQ $Merged_Rv_FastQ"

    MY_BWA_PREFIX=$HOME/results/$SLX.$SPECIES.v1/BWA/$Barcode/$SLX.$Barcode.$Cell.1K.bwa
    echo -e "bwa mem -M $TRIMMED_FASTQ_FILES | samblaster | samtools \n"
    time bwa mem \
        -M \
        -t $NT \
        -R "@RG\tID:$SLX\tPL:illumina\tPU:run\tLB:$SLX\tSM:$Barcode\tCN:CamObsGynae" \
        $BWA_INDEX_BASE.fa \
        $TRIMMED_FASTQ_FILES 2> $MY_BWA_PREFIX.log | \
        time $HOME/Install/samblaster/samblaster -M 2> $MY_BWA_PREFIX.samblaster.log | \
        time samtools view -Sb - > $MY_BWA_PREFIX.samblaster.bam
