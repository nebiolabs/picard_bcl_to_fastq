#!/bin/sh

#argument: path to the samplesheet in a run folder
PICARD_PATH=/mnt/ngswork/galaxy/sw/picard-tools-1.94/
READ_STRUCTURE=36T8B
CPU_COUNT=32
JAVA_OPTS=-Xmx50g

if [ $# -eq 0 ]; then
	echo "specify the sample sheet"
	exit 1
fi

sample_sheet=`readlink -f "${1}"`
echo "Sample sheet: ${sample_sheet}"
run_path=`dirname "${sample_sheet}"`
echo "Run path: ${run_path}"
run_barcode=`echo "${run_path}" | rev | cut -f1 -d'/' | rev | cut -f2 -d'-'`
echo "Barcode: ${run_barcode}"

if [ ! -d ${run_path}/fastq ] ; then
	mkdir ${run_path}/fastq
fi

echo 'barcode_sequence_1	barcode_name	library_name' > barcode_params.txt
perl -nle 'print "$3\t$2\t$1" if /^([^,]+)(?:[^,]*,){4}([^,]+),([GCAT]+),.*$/' "${sample_sheet}" >> barcode_params.txt
echo 'OUTPUT_PREFIX	BARCODE_1' > multiplex_params.txt
echo 'unassigned	N' >> multiplex_params.txt
perl -nle 'print "$1\t$3" if /^([^,]+)(?:[^,]*,){4}([^,]+),([GCAT]+),.*$/' "${sample_sheet}" >> multiplex_params.txt

cd "${run_path}/fastq"
java  $JAVA_OPTS -jar $PICARD_PATH/ExtractIlluminaBarcodes.jar  NUM_PROCESSORS=$CPU_COUNT READ_STRUCTURE=$READ_STRUCTURE LANE=001 BASECALLS_DIR=${run_path}/Data/Intensities/BaseCalls METRICS_FILE=barcode_metrics.txt BARCODE_FILE=${run_path}/barcode_params.txt
java  $JAVA_OPTS -jar $PICARD_PATH/IlluminaBasecallsToFastq.jar NUM_PROCESSORS=$CPU_COUNT READ_STRUCTURE=$READ_STRUCTURE RUN_BARCODE=$run_barcode LANE=1 BASECALLS_DIR=${run_path}/Data/Intensities/BaseCalls MULTIPLEX_PARAMS=${run_path}/multiplex_params.txt
