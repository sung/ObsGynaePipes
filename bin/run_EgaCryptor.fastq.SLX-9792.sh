SLX="SLX-9792"
BARCODES="D701_D502
D701_D504
D701_D506
D701_D508
D709_D502
D709_D504
D709_D506
D709_D508
D710_D502
D710_D504
D710_D506
D710_D508
D711_D508
D712_D502
D712_D506
"

for i in $BARCODES; do echo -e "SLX:$SLX Barcode:$i"; java -jar ~/Install/EgaCryptor/EgaCryptor.jar -file /remote/NAS1/data/fastq/$SLX/$SLX.$i*.fq.gz; done
