#!/bin/bash

cufflinks -p 8 --library-type fr-firststrand --GTF-guide /whale-data/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf --output-dir /whale-data/ssg29/RNA-Seq/Results/Pilot1_v6/Cufflink/R4 /whale-data/ssg29/RNA-Seq/Results/Pilot1_v6/TopHat/R4/accepted_hits.bam 

