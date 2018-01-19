#!/bin/bash

SLX="SLX-9168"
BARCODES="D701_D502 D701_D504 D701_D506 D701_D508 D702_D502 D702_D504 D702_D506 D702_D508 D707_D502 D707_D503 D707_D506 D707_D508 D708_D502 D708_D504 D708_D508 D709_D502 D709_D504 D709_D506 D709_D508 D710_D502 D710_D504 D710_D506 D710_D508 D711_D502 D711_D504 D711_D506 D711_D508"
for i in $BARCODES; do for j in `ls /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz`; do FILE=`echo $j | cut -d'/' -f 7`; k=`cat $j.md5`; l=`cat $j.gpg.md5`; printf "$SLX\t$i\t/RNA-Seq/fastq/$FILE.gpg\t$l\t$k\n"; done; done;

SLX="SLX-9169"
BARCODES="D704_D502 D704_D504 D704_D506 D704_D508 D705_D502 D705_D504 D705_D506 D705_D508 D706_D502 D706_D503 D706_D505 D706_D506 D706_D508 D708_D502 D709_D502 D709_D504 D709_D506 D709_D508 D710_D502 D710_D504 D710_D506 D710_D508 D711_D502 D711_D504 D711_D506 D711_D508 D712_D502 D712_D504 D712_D506 D712_D508"
for i in $BARCODES; do for j in `ls /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz`; do FILE=`echo $j | cut -d'/' -f 7`; k=`cat $j.md5`; l=`cat $j.gpg.md5`; printf "$SLX\t$i\t/RNA-Seq/fastq/$FILE.gpg\t$l\t$k\n"; done; done;

SLX="SLX-10281"
BARCODES="D701_D502 D701_D504 D701_D506 D701_D508 D702_D502 D702_D507 D702_D508 D711_D502 D711_D504 D711_D506 D712_D502 D712_D504 D712_D506 D712_D508"
for i in $BARCODES; do for j in `ls /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz`; do FILE=`echo $j | cut -d'/' -f 7`; k=`cat $j.md5`; l=`cat $j.gpg.md5`; printf "$SLX\t$i\t/RNA-Seq/fastq/$FILE.gpg\t$l\t$k\n"; done; done;

SLX="SLX-10402"
BARCODES="D704_D502 D704_D504 D704_D506 D704_D508 D705_D502 D705_D504 D705_D506 D705_D508 D706_D502 D706_D504 D706_D506 D706_D508 D707_D502 D707_D504 D707_D506 D707_D508"
for i in $BARCODES; do for j in `ls /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz`; do FILE=`echo $j | cut -d'/' -f 7`; k=`cat $j.md5`; l=`cat $j.gpg.md5`; printf "$SLX\t$i\t/RNA-Seq/fastq/$FILE.gpg\t$l\t$k\n"; done; done;

SLX="SLX-9792"
BARCODES="D701_D502 D701_D504 D701_D506 D701_D508 D709_D502 D709_D504 D709_D506 D709_D508 D710_D502 D710_D504 D710_D506 D710_D508 D711_D508 D712_D502 D712_D506"
for i in $BARCODES; do for j in `ls /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz`; do FILE=`echo $j | cut -d'/' -f 7`; k=`cat $j.md5`; l=`cat $j.gpg.md5`; printf "$SLX\t$i\t/RNA-Seq/fastq/$FILE.gpg\t$l\t$k\n"; done; done;

SLX="SLX-10284"
BARCODES="D701_D502 D701_D504 D701_D506 D701_D508 D702_D502 D702_D504 D702_D506 D702_D508 D704_D502 D704_D504 D704_D508"
for i in $BARCODES; do for j in `ls /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz`; do FILE=`echo $j | cut -d'/' -f 7`; k=`cat $j.md5`; l=`cat $j.gpg.md5`; printf "$SLX\t$i\t/RNA-Seq/fastq/$FILE.gpg\t$l\t$k\n"; done; done;

SLX="SLX-10285"
BARCODES="D705_D502 D705_D504 D705_D506 D705_D508 D706_D502 D706_D504 D706_D506 D706_D508 D707_D502 D707_D506 D707_D508"
for i in $BARCODES; do for j in `ls /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz`; do FILE=`echo $j | cut -d'/' -f 7`; k=`cat $j.md5`; l=`cat $j.gpg.md5`; printf "$SLX\t$i\t/RNA-Seq/fastq/$FILE.gpg\t$l\t$k\n"; done; done;

SLX="SLX-10283"
BARCODES="D709_D502 D709_D504 D709_D506 D709_D508 D710_D502 D710_D504 D710_D508 D712_D506 D712_D508"
for i in $BARCODES; do for j in `ls /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz`; do FILE=`echo $j | cut -d'/' -f 7`; k=`cat $j.md5`; l=`cat $j.gpg.md5`; printf "$SLX\t$i\t/RNA-Seq/fastq/$FILE.gpg\t$l\t$k\n"; done; done;

SLX="SLX-10287"
BARCODES="D701_D502 D701_D504 D701_D506 D701_D508 D702_D502 D702_D504 D702_D506 D702_D508 D704_D502 D704_D508 D705_D502 D705_D504 D705_D506 D706_D502"
for i in $BARCODES; do for j in `ls /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz`; do FILE=`echo $j | cut -d'/' -f 7`; k=`cat $j.md5`; l=`cat $j.gpg.md5`; printf "$SLX\t$i\t/RNA-Seq/fastq/$FILE.gpg\t$l\t$k\n"; done; done;
