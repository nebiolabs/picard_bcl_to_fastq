#!/bin/sh

#argument: path to the samplesheet in a run folder
PICARD_PATH=/mnt/ngswork/galaxy/sw/picard-tools-1.97/
CPU_COUNT=`grep -i processor /proc/cpuinfo | wc -l`
echo "Running on ${CPU_COUNT} cpus"
ulimit -n 4096

JAVA_OPTS=-Xmx50g 

if [ $# -eq 0 ]; then
	echo "specify the sample sheet"
	exit 1
fi

sample_sheet=`readlink -f "${1}"`
echo "Sample sheet: ${sample_sheet}"
run_path=`dirname "${sample_sheet}"`
echo "Run path: ${run_path}"

barcode_params="${run_path}/barcode_params.txt"
multiplex_params="${run_path}/multiplex_params.txt" 

run_barcode=`echo "${run_path}" | rev | cut -f1 -d'/' | rev | cut -f2 -d'-'`
echo "Barcode: ${run_barcode}"

if [ ! -d "${run_path}/fastq" ] ; then
	mkdir "${run_path}/fastq"
fi

echo 'barcode_sequence_1	barcode_name	library_name' > "${run_path}/barcode_params.txt"
perl -nle 'print "$3\t$2\t$1" if /^([^,]+)(?:[^,]*,){3,4}([^,]+),([GCAT]+),.*$/' "${sample_sheet}" >> "${run_path}/barcode_params.txt"
echo 'OUTPUT_PREFIX	BARCODE_1' > ${multiplex_params}
echo 'unassigned	N' >> ${multiplex_params}
perl -nle 'print "$1\t$3" if /^([^,]+)(?:[^,]*,){3,4}([^,]+),([GCAT]+),.*$/' "${sample_sheet}" >> ${multiplex_params}

read_cycles=`perl -nle 'print "$1" if /^(\d+),*\s*$/' "${sample_sheet}"`
read1_cycles=`echo ${read_cycles} | cut -d' ' -f1`
read2_cycles=`echo ${read_cycles} | cut -d' ' -f2 -s`

bc_cycles=`perl -nle 'if (/^([^,]+)(?:[^,]*,){3,4}([^,]+),([GCAT]+),.*$/) { print length($3); exit}' "${sample_sheet}"`
READ_STRUCTURE="${read1_cycles}T${bc_cycles}B"
if [ -n "${read2_cycles}" ] ; then
	READ_STRUCTURE="${READ_STRUCTURE}${read2_cycles}T"
fi 
echo "read structure: ${READ_STRUCTURE}"

barcode_count=$((`wc -l "${barcode_params}" | cut -d ' ' -f1` - 1))

if [ $barcode_count -eq 0 ]; then
	echo "Failed to find any barcodes in ${sample_sheet}" 1>&2
	exit 1
fi

cd "${run_path}/fastq"

java  $JAVA_OPTS -jar $PICARD_PATH/ExtractIlluminaBarcodes.jar  MAX_NO_CALLS=2 MIN_MISMATCH_DELTA=2 MAX_MISMATCHES=1 NUM_PROCESSORS=$CPU_COUNT READ_STRUCTURE=$READ_STRUCTURE LANE=001 BASECALLS_DIR="${run_path}/Data/Intensities/BaseCalls" METRICS_FILE=barcode_metrics.txt BARCODE_FILE="${barcode_params}"
java  $JAVA_OPTS -jar $PICARD_PATH/IlluminaBasecallsToFastq.jar NUM_PROCESSORS=$CPU_COUNT READ_STRUCTURE=$READ_STRUCTURE RUN_BARCODE=$run_barcode LANE=1 BASECALLS_DIR="${run_path}/Data/Intensities/BaseCalls" MULTIPLEX_PARAMS="${multiplex_params}"
