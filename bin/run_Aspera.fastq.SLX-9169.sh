SLX="SLX-9169"
BARCODES="D704_D502
D704_D504
D704_D506
D704_D508
D705_D502
D705_D504
D705_D506
D705_D508
D706_D502
D706_D503
D706_D505
D706_D506
D706_D508
D708_D502
D709_D502
D709_D504
D709_D506
D709_D508
D710_D502
D710_D504
D710_D506
D710_D508
D711_D502
D711_D504
D711_D506
D711_D508
D712_D502
D712_D504
D712_D506
D712_D508
"
export ASPERA_SCP_PASS=CZ9BfwsP


for i in $BARCODES; do echo -e "SLX:$SLX Barcode:$i"; time ascp -QT -l300M -L- /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz.gpg* ega-box-651@fasp.ega.ebi.ac.uk:/RNA-Seq/fastq/; done
