####################
## WGBS on-target ##
####################
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/NKX1-2/NKX1-2.up.1500.down.500.bed > ${i%.bed}.NKX1-2.up.1500.down.500.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/MAB21L1/MAB21L1.top2.grch37.bed > ${i%.bed}.MAB21L1.top2.500bp.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/MAB21L1/MAB21L1.top2.grch37.enstTiling500.bed > ${i%.bed}.MAB21L1.top2.enstTiling.500bp.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/Chr7.DMR/chr7.top2.bed > ${i%.bed}.chr7.DMR.top2.500bp.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/VIPR2/VIPR2.top5k.grch37.bed > ${i%.bed}.VIPR2.DMR.5kb.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/C1orf216/C1orf216.up.1500.down.500.grch37.bed > ${i%.bed}.C1orf216.up.1500.down.500.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/TICAM2/TICAM2.grch37.bed > ${i%.bed}.TICAM2.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr.bed > ${i%.bed}.CSMD1.DMR.225kb.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/CSMD1/CSMD1.placenta.putative.promoter.4K.GRCh37.bed > ${i%.bed}.CSMD1.placenta.putative.promoter.4K.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/CSMD1/CSMD1.legacy.promo.4K.GRCh37.bed > ${i%.bed}.CSMD1.legacy.promoter.4K.bed ; done
#for i in `grep bed  /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv | cut -d, -f1`; do intersectBed -header -a $i -b ~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr+2exons.bed > ${i%.bed}.CSMD1.DMR.247kb.bed ; done
for i in `awk -F',' 'NR>1{print $1}' /scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv`; do intersectBed -header -a $i -b ~/Pipelines/data/SMS/SMS.GRCh37.1kb.up.down.UCSC.bed > ${i%.bed.gz}.SMS.1kb.up.down.bed ; done

##########################
## SureSelect on-target ##
##########################
#for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/NKX1-2/NKX1-2.up.1500.down.500.bed > ${i%.bed}.NKX1-2.up.1500.down.500.bed ; done
#for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/MAB21L1/MAB21L1.top2.grch37.bed > ${i%.bed}.MAB21L1.top2.500bp.bed ; done 
#for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/Chr7.DMR/chr7.top2.bed > ${i%.bed}.chr7.DMR.top2.500bp.bed ; done
#for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/VIPR2/VIPR2.top5k.grch37.bed > ${i%.bed}.VIPR2.DMR.5kb.bed ; done
#for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/C1orf216/C1orf216.up.1500.down.500.grch37.bed > ${i%.bed}.C1orf216.up.1500.down.500.bed ; done
#for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/TICAM2/TICAM2.grch37.bed > ${i%.bed}.TICAM2.bed ; done
#for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr.bed > ${i%.bed}.CSMD1.DMR.225kb.bed ; done
#for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/CSMD1/CSMD1.legacy.promo.4K.GRCh37.bed > ${i%.bed}.CSMD1.legacy.promoter.4K.bed ; done
#for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr+2exons.bed > ${i%.bed}.CSMD1.DMR.247kb.bed ; done
for i in `ls /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.CG.bed`; do intersectBed -header -a $i -b ~/Pipelines/data/SMS/SMS.GRCh37.1kb.up.down.UCSC.bed > ${i%.bed}.SMS.1kb.up.down.bed ; done

###########################################################
### Concordance of CpG calls between WGBS and SureSelect ##
###########################################################
# NKX1-2 (chr10)
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.CG.NKX1-2.up.1500.down.500.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.CG.NKX1-2.up.1500.down.500.bed -wa -wb -s > /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.CG.NKX1-2.up.1500.down.500.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.CG.NKX1-2.up.1500.down.500.bed -wa -wb -s > /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.CG.NKX1-2.up.1500.down.500.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A002/SLX-8075.SLX-8077.SLX-8081.A002.cpg.filtered.CG.NKX1-2.up.1500.down.500.bed -wa -wb -s > /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.CG.NKX1-2.up.1500.down.500.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A004/SLX-8075.SLX-8077.SLX-8081.A004.cpg.filtered.CG.NKX1-2.up.1500.down.500.bed -wa -wb -s > /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.top.concordance.with.wgbs.bed
# MAB21L1 (chr13)
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.CG.MAB21L1.top2.500bp.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.CG.MAB21L1.top2.500bp.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.CG.MAB21L1.top2.500bp.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.CG.MAB21L1.top2.500bp.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.CG.MAB21L1.top2.500bp.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A002/SLX-8075.SLX-8077.SLX-8081.A002.cpg.filtered.CG.MAB21L1.top2.500bp.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.CG.MAB21L1.top2.500bp.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A004/SLX-8075.SLX-8077.SLX-8081.A004.cpg.filtered.CG.MAB21L1.top2.500bp.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.top.concordance.with.wgbs.bed
# chr7.DMR (chr7)
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.CG.chr7.DMR.top2.500bp.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.CG.chr7.DMR.top2.500bp.bed -wa -wb -s > /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.CG.chr7.DMR.top2.500bp.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.CG.chr7.DMR.top2.500bp.bed -wa -wb -s > /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.CG.chr7.DMR.top2.500bp.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A002/SLX-8075.SLX-8077.SLX-8081.A002.cpg.filtered.CG.chr7.DMR.top2.500bp.bed -wa -wb -s > /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.CG.chr7.DMR.top2.500bp.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A004/SLX-8075.SLX-8077.SLX-8081.A004.cpg.filtered.CG.chr7.DMR.top2.500bp.bed -wa -wb -s > /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.top.concordance.with.wgbs.bed
# CSMD1 (chr8)
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.CG.CSMD1.DMR.225kb.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.CG.CSMD1.DMR.225kb.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.CG.CSMD1.DMR.225kb.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.CG.CSMD1.DMR.225kb.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.CG.CSMD1.DMR.225kb.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A002/SLX-8075.SLX-8077.SLX-8081.A002.cpg.filtered.CG.CSMD1.DMR.225kb.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.CG.CSMD1.DMR.225kb.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A004/SLX-8075.SLX-8077.SLX-8081.A004.cpg.filtered.CG.CSMD1.DMR.225kb.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.top.concordance.with.wgbs.bed
# C1orf216 (chr1)
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.CG.C1orf216.up.1500.down.500.grch37.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.CG.C1orf216.up.1500.down.500.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.CG.C1orf216.up.1500.down.500.grch37.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.CG.C1orf216.up.1500.down.500.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.CG.C1orf216.up.1500.down.500.grch37.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A002/SLX-8075.SLX-8077.SLX-8081.A002.cpg.filtered.CG.C1orf216.up.1500.down.500.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.CG.C1orf216.up.1500.down.500.grch37.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A004/SLX-8075.SLX-8077.SLX-8081.A004.cpg.filtered.CG.C1orf216.up.1500.down.500.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.top.concordance.with.wgbs.bed
#TICAM2 (chr5)
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.CG.TICAM2.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.CG.TICAM2.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL03/SLX-10409.HAL03.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.CG.TICAM2.bed -b /scratch/obsgynae/POPS/results/SLX-8074.SLX-8080.Homo_sapiens.v1/BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.CG.TICAM2.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL46/SLX-10409.HAL46.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.CG.TICAM2.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A002/SLX-8075.SLX-8077.SLX-8081.A002.cpg.filtered.CG.TICAM2.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL01/SLX-10409.HAL01.cpg.filtered.top.concordance.with.wgbs.bed
#bedtools intersect -a /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.CG.TICAM2.bed -b /scratch/obsgynae/POPS/results/SLX-8075.SLX-8077.SLX-8081.Homo_sapiens.v1/BisSNP/A004/SLX-8075.SLX-8077.SLX-8081.A004.cpg.filtered.CG.TICAM2.bed -wa -wb -s >> /scratch/obsgynae/POPS/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL33/SLX-10409.HAL33.cpg.filtered.top.concordance.with.wgbs.bed