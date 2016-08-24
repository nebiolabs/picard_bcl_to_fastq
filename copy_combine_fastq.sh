#!/bin/bash

flowcell=`/opt/bio/bin/xmllint --xpath " //Flowcell/text()" RunInfo.xml`
run_date=`/opt/bio/bin/xmllint --xpath " //Date/text()" RunInfo.xml`

echo "Copying data to combined_fastq"
mkdir -p combined_fastq

for f in `cut -f3  lane1_barcode_params.txt  | tail -n +2 | tr ' ' '-' | tr "\n" " "` unassigned ;
do
  for r in 1 2 ;
  do
    echo "combining " *_$f.$r.fastq
    cat fastq/*/*_$f.$r.fastq | /mnt/galaxy/data/galaxy/sw/bin/pigz -p 10 > "combined_fastq/$f.$r.fastq.gz";
  done;
done;