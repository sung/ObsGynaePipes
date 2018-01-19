#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 27/Aug/2015
# Last modified: 27/Aug/2015
# Optimised and customised to run at the Darwin HPC

source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

MIRBASE_VERSION=20
DEST=$HOME/data/Annotation/miRbase.$MIRBASE_VERSION
mkdir_unless $DEST
chdir $DEST

time wget --directory-prefix=$DEST ftp://mirbase.org/pub/mirbase/$MIRBASE_VERSION/aliases.txt.gz
time wget --directory-prefix=$DEST ftp://mirbase.org/pub/mirbase/$MIRBASE_VERSION/hairpin.fa.gz
time wget --directory-prefix=$DEST ftp://mirbase.org/pub/mirbase/$MIRBASE_VERSION/mature.fa.gz
time wget --directory-prefix=$DEST ftp://mirbase.org/pub/mirbase/$MIRBASE_VERSION/miFam.dat.gz
time wget --directory-prefix=$DEST ftp://mirbase.org/pub/mirbase/$MIRBASE_VERSION/miRNA.dat.gz
time wget --directory-prefix=$DEST ftp://mirbase.org/pub/mirbase/$MIRBASE_VERSION/miRNA.str.gz
time wget --directory-prefix=$DEST ftp://mirbase.org/pub/mirbase/$MIRBASE_VERSION/miRNA.xls.gz
time wget --directory-prefix=$DEST ftp://mirbase.org/pub/mirbase/$MIRBASE_VERSION/organisms.txt.gz
time wget --directory-prefix=$DEST ftp://mirbase.org/pub/mirbase/$MIRBASE_VERSION/genomes/hsa.gff3
# above based on GRCh38

for i in `ls $DEST/*.gz`;do  echo $i; gunzip -f $i; done

# hsa
awk 'BEGIN{ ORS = ""; RS = ">"; FS="\n" }$1~"hsa"{ print ">" $0 }' $DEST/hairpin.fa | awk '{print $1}' > $DEST/hsa.hairpin.fa
awk 'BEGIN{ ORS = ""; RS = ">"; FS="\n" }$1~"hsa"{ print ">" $0 }' $DEST/mature.fa | awk '{print $1}' > $DEST/hsa.mature.fa
for i in `grep Primates $DEST/organisms.txt | grep -v hsa | cut -f1`;do awk "BEGIN{ORS = \"\"; RS = \">\"; FS=\"\n\"}\$1~\"$i\"{ print \">\" \$0 }" $DEST/mature.fa; done | awk '{print $1}' > $DEST/hsa.mature.homologues.fa

# mmu
awk 'BEGIN{ ORS = ""; RS = ">"; FS="\n" }$1~"mmu"{ print ">" $0 }' $DEST/hairpin.fa | awk '{print $1}'> $DEST/mmu.hairpin.fa
awk 'BEGIN{ ORS = ""; RS = ">"; FS="\n" }$1~"mmu"{ print ">" $0 }' $DEST/mature.fa | awk '{print $1}'> $DEST/mmu.mature.fa
for i in `grep Rodentia $DEST/organisms.txt | grep -v mmu | cut -f1`;do awk "BEGIN{ORS = \"\"; RS = \">\"; FS=\"\n\"}\$1~\"$i\"{ print \">\" \$0 }" $DEST/mature.fa; done | awk '{print $1}' > $DEST/mmu.mature.homologues.fa
