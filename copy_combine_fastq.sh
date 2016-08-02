#!/bin/bash
if [[ -z  "$1" ]]; then
        echo "please specify a user folder to import (e.g. langhorst@neb.com)"
        exit 1
fi
echo "user: $1"

galaxy_user_import_base='/mnt/galaxy/tmp/users'
flowcell=`/opt/bio/bin/xmllint --xpath " //Flowcell/text()" RunInfo.xml`
run_date=`/opt/bio/bin/xmllint --xpath " //Date/text()" RunInfo.xml`
dest_dir="${galaxy_user_import_base}/$1/${run_date}_${flowcell}"

if  [ ! -d  "$galaxy_user_import_base/$1" ]; then
        echo "Could not find $1 in $galaxy_user_import_base"
        exit 1
fi
echo "Copying data to $dest_dir"
mkdir -p $dest_dir

for f in `cut -f3  lane1_barcode_params.txt  | tail -n +2 | tr ' ' '-' | tr "\n" " "` unassigned ;
do
  for r in 1 2 ;
  do
    echo "combining " *_$f.$r.fastq
    cat fastq/*/*_$f.$r.fastq | /mnt/galaxy/data/galaxy/sw/bin/pigz -p 28 > "$dest_dir/$f.$r.fastq.gz";
  done;
done;