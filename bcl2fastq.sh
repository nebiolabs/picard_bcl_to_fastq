#!/bin/sh

#argument: path to the samplesheet in a run folder

#Path setup
PICARD_PATH=/mnt/ngswork/galaxy/sw/picard-tools-1.97/
CPU_COUNT=`grep -i processor /proc/cpuinfo | wc -l`

#demultiplexer settings
MAX_MISMATCHES=0
MAX_NO_CALLS=0
MIN_MISMATCH_DELTA=2

#Java settings
JAVA_OPTS=-Xmx100g 

echo "Running on ${CPU_COUNT} cpus"

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

if [ ! -d "${run_path}/fastq" ] ; then
	mkdir "${run_path}/fastq"
fi

if [ ! -x `which xmllint` ] ; then
	echo "xmllint program is not found on the path, cannot continue."
	exit 1
fi

read1_cycles=`echo 'cat //Read[@Number="1"]/@NumCycles' | xmllint -shell "${run_path}/RunInfo.xml"  | sed -n 3p | sed s/.*=// | sed s/\"//g`
bc_cycles=`echo 'cat //Read[@Number="2"]/@NumCycles' | xmllint -shell "${run_path}/RunInfo.xml"  | sed -n 3p | sed s/.*=// | sed s/\"//g`
read2_cycles=`echo 'cat //Read[@Number="3"]/@NumCycles' | xmllint -shell "${run_path}/RunInfo.xml"  | sed -n 3p | sed s/.*=// | sed s/\"//g`

read_structure="${read1_cycles}T"

if [ "$bc_cycles"+0 -gt "0" ]; then
	read_structure="${read_structure}${bc_cycles}B"
fi

if [ -n "${read2_cycles}" ] ; then
	read_structure="${read_structure}${read2_cycles}T"
fi 
echo "read structure: ${read_structure}"

cd "${run_path}/fastq"

metrics_name="${MAX_NO_CALLS}nc_${MIN_MISMATCH_DELTA}mmd_${MAX_MISMATCHES}mis_bc_metrics.txt"

lanecount=`echo 'cat //FlowcellLayout/@LaneCount' | xmllint -shell "${run_path}/RunInfo.xml"  | sed -n 3p | sed s/.*=// | sed s/\"//g`
if [ "$lanecount" -eq "1" ]; then
	echo "Detected miseq format"
	is_miseq=true
else
	is_miseq=false
fi

for i in `seq 1 ${lanecount}`
do
	echo processing lane ${i} 

	multiplex_params="${run_path}/lane${i}_multiplex_params.txt" 
	barcode_params="${run_path}/lane${i}_barcode_params.txt"

	echo 'barcode_sequence_1	barcode_name	library_name' > "${barcode_params}"
	regex=
	if $is_miseq ; then
		regex="/^([^,]+)(?:[^,]*,){3,4}([^,]+),([GCAT]+),.*$/"
		perl -nle "print \"\$3\t\$2\t\$1\" if ${regex}" "${sample_sheet}" >> "${barcode_params}"
	else
		regex="/^(?:[^,]+),${i},([^,]+),(?:[^,]*,)([GCAT]*),(?:[^,]*,){4}\w+\s*$/"
		perl -nle "print ((\$2 || 'N').\"\t\$1\t\$1\") if ${regex}" "${sample_sheet}" >> "${barcode_params}"
	fi

	barcode_count=$((`wc -l "${barcode_params}" | cut -d ' ' -f1` - 1))
	if [ $barcode_count -eq 0 ]; then
		echo "Warning: Failed to find any barcodes in ${sample_sheet}" 1>&2
	fi
		
	echo 'OUTPUT_PREFIX	BARCODE_1' > ${multiplex_params}
	if [ "${barcode_count}" -gt "1" ] || [ $barcode_count -eq 0 ]; then
		echo "L${i}_unassigned	N" >> ${multiplex_params}
	fi

	if $is_miseq ; then
		perl -nle "print \"\$1\t\$3\" if ${regex}" "${sample_sheet}" >> "${multiplex_params}"
	else
		perl -nle "print \"L${i}_\$1\t\$2\" if ${regex} " "${sample_sheet}" >> "${multiplex_params}"
	fi
	#cat $barcode_params
	#cat $multiplex_params
	
	if [ $barcode_count -gt 0 ]; then

  	    java  $JAVA_OPTS -jar $PICARD_PATH/ExtractIlluminaBarcodes.jar \
		MAX_NO_CALLS=$MAX_NO_CALLS MIN_MISMATCH_DELTA=$MIN_MISMATCH_DELTA \
		MAX_MISMATCHES=$MAX_MISMATCHES NUM_PROCESSORS=$CPU_COUNT \
		read_structure=$read_structure \
		LANE=${i} \
		BASECALLS_DIR="${run_path}/Data/Intensities/BaseCalls" \
		METRICS_FILE="L${i}_${metrics_name}" BARCODE_FILE="${barcode_params}"
	fi

	java  $JAVA_OPTS -jar $PICARD_PATH/IlluminaBasecallsToFastq.jar \
		NUM_PROCESSORS=$CPU_COUNT read_structure=$read_structure \
		RUN_BARCODE=$run_barcode LANE=${i} \
		BASECALLS_DIR="${run_path}/Data/Intensities/BaseCalls" MULTIPLEX_PARAMS="${multiplex_params}"
done
