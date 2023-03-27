
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

run_path = workflow.launchDir.toString()
path_to_java = params.path_to_java
path_to_picard = params.path_to_picard

process writeParamsFiles {
    cpus 1

    input:
        tuple val(lane), val(sample_sheet_data)

    output:
        tuple val(lane), file("*_library_params*.txt"), file("*_barcode_params*.txt"), emit: individual_params_files

    shell:
    // below is the format of barcode_params
    // barcode_sequence_1      barcode_sequence_2      barcode_name    library_name
    // ACTTCTGC        AGTCCCGG        ACTTCTGCAGTCCCGG        FS_500ng_10minFS_lot0923_Loop_KapaHiFi_x4_NoSPRI_E6447_A1

    // below is the format of library_params
    // BARCODE_1       BARCODE_2       OUTPUT  LIBRARY_NAME    SAMPLE_ALIAS
    // ACTTCTGC        AGTCCCGG        L1_FS_500ng_10minFS_lot0923_Loop_KapaHiFi_x4_NoSPRI_E6447_A1.bam        FS_500ng_10minFS_lot0923_Loop_KapaHiFi_x4_NoSPRI_E6447_A1       FS_500ng_10minFS_lot0923_Loop_KapaHiFi_x

    // need to account for single or dual index

    '''
    echo -e "!{sample_sheet_data.index}\\t!{sample_sheet_data.index2}\\t!{sample_sheet_data.index}!{sample_sheet_data.index2}\\t!{sample_sheet_data.Sample_ID}" | sed 's/\\t\\t/\\t/g' > !{sample_sheet_data.Sample_ID}_barcode_params.!{lane}.txt
    echo -e "!{sample_sheet_data.index}\\t!{sample_sheet_data.index2}\\tL!{lane}_!{sample_sheet_data.Sample_ID}.bam\\t!{sample_sheet_data.Sample_ID}\\t!{sample_sheet_data.Sample_ID}" | sed 's/\\t\\t/\\t/g' > !{sample_sheet_data.Sample_ID}_library_params.!{lane}.txt
    '''
}

process combineParamsFiles{
    cpus 1
    publishDir "${run_path}", mode: 'copy', overwrite: true

    input:
        tuple val(lane), file(individual_library_params), file(individual_barcode_params)

    output:
        tuple val(lane), file("barcode_params*.txt"), emit: barcode_params_file
        tuple val(lane), file("library_params*.txt"), emit: library_params_file

    shell:
    // below is the format of barcode_params
    // barcode_sequence_1      barcode_sequence_2      barcode_name    library_name
    // ACTTCTGC        AGTCCCGG        ACTTCTGCAGTCCCGG        FS_500ng_10minFS_lot0923_Loop_KapaHiFi_x4_NoSPRI_E6447_A1

    // below is the format of library_params
    // BARCODE_1       BARCODE_2       OUTPUT  LIBRARY_NAME    SAMPLE_ALIAS
    // ACTTCTGC        AGTCCCGG        L1_FS_500ng_10minFS_lot0923_Loop_KapaHiFi_x4_NoSPRI_E6447_A1.bam        FS_500ng_10minFS_lot0923_Loop_KapaHiFi_x4_NoSPRI_E6447_A1       FS_500ng_10minFS_lot0923_Loop_KapaHiFi_x

    '''
    if [[ $(awk '{print NF}' !{individual_barcode_params[0]} | uniq) -eq 3 ]]; then # single index
        cat <(echo -e "barcode_sequence_1\\tbarcode_name\\tlibrary_name") *_barcode_params.!{lane}.txt > barcode_params.!{lane}.txt
        cat <(echo -e "BARCODE_1\\tOUTPUT\\tLIBRARY_NAME\\tSAMPLE_ALIAS") *_library_params.!{lane}.txt > library_params.!{lane}.txt
        cat <(echo -e "N\\tL!{lane}_unmatched.bam\\tunmatched\\tunmatched") >> library_params.!{lane}.txt
    else # dual index
        cat <(echo -e "barcode_sequence_1\\tbarcode_sequence_2\\tbarcode_name\\tlibrary_name") *_barcode_params.!{lane}.txt > barcode_params.!{lane}.txt
        cat <(echo -e "BARCODE_1\\tBARCODE_2\\tOUTPUT\\tLIBRARY_NAME\\tSAMPLE_ALIAS") *_library_params.!{lane}.txt > library_params.!{lane}.txt
        cat <(echo -e "N\\tN\\tL!{lane}_unmatched.bam\\tunmatched\\tunmatched") >> library_params.!{lane}.txt
    fi
    '''

}

process extractIlluminaBarcodes{
    cpus 16

    conda "picard=2.27.4"
    publishDir "${run_path}/tile_bams", mode: 'symlink', overwrite: true

    input:
        tuple val(lane), file(barcode_params)

    output:
        tuple val(lane), file('L*.txt'), emit: demux_stats_file

    shell:
    ''' 
    !{path_to_java}/java -Xmx65536m -jar !{path_to_picard}/picard.jar ExtractIlluminaBarcodes \
            MAX_NO_CALLS=!{params.max_no_calls} MIN_MISMATCH_DELTA=!{params.min_mismatch_delta} \
            MAX_MISMATCHES=!{params.max_mismatches} NUM_PROCESSORS=!{task.cpus} \
            read_structure=!{params.read_structure} \
            LANE=!{lane} \
            BASECALLS_DIR=\"!{run_path}/Data/Intensities/BaseCalls\" \
            METRICS_FILE=\"L!{lane}_!{params.max_no_calls}nc_!{params.min_mismatch_delta}mmd_!{params.max_mismatches}mis_bc_metrics.txt\" \
            BARCODE_FILE=\"!{barcode_params}\" \
            COMPRESS_OUTPUTS=true \
            DISTANCE_MODE=LENIENT_HAMMING
    '''
}

process updateReadCounts {
    cpus 1

    input:
        val(_file_checked)

    shell:
    '''
    ruby "!{workflow.projectDir}/post_demuxed_reads.rb" "!{run_path}" "!{params.pool_id}" "!{params.read_structure.count("T")}"
    '''
}

process illuminaBasecallsToSam {
    cpus 4

    conda "picard=2.27.4"
    publishDir "${run_path}/tile_bams/L_${lane}_${tile}", mode: 'symlink', overwrite: true

    input:
        tuple val(lane), file(library_params_file), val(tile), file(demux_metrics)

    output:
        path("*.bam"), emit: tile_bams // file("*.bam") doesn't work with emit:
        val("${run_path}/tile_bams/${lane}_${tile}"), emit: tile_dir

    shell:
    '''
    !{path_to_java}/java -Xmx16384m -jar !{path_to_picard}/picard.jar IlluminaBasecallsToSam \
                LIBRARY_PARAMS=\"!{library_params_file}\" \
                BASECALLS_DIR=\"!{run_path}/Data/Intensities/BaseCalls\" \
                NUM_PROCESSORS=!{task.cpus} \
                READ_STRUCTURE=!{params.read_structure} \
                RUN_BARCODE=!{params.machine}:!{params.machine}:!{params.flowcell} \
                PROCESS_SINGLE_TILE=!{tile} \
                LANE=!{lane} \
                MAX_READS_IN_RAM_PER_TILE=3000000 \
                MAX_RECORDS_IN_RAM=3000000 \
                PLATFORM=ILLUMINA \
                INCLUDE_BC_IN_RG_TAG=true \
                TMP_DIR=/state/partition1/sge_tmp \
                BARCODE_POPULATION_STRATEGY=ALWAYS
    '''
}

process findTiles {
    cpus 1

    input:
        tuple val(lane), file(tile_info_file)

    output:
        tuple val(lane), stdout, emit: tiles

    shell:
    '''
    echo 'cat //Tile/text()' |
    xmllint -shell !{tile_info_file} |
    grep -P '^\\d' |
    sort -n |
    awk '/^[0-9]+_[0-9]+$/ {print $1} !/_/ {print "1_"$1}' |
    grep !{lane}_ |
    cut -f 2 -d "_"
    '''
}

process checkReadCounts {
    cpus 1
    cache false
    tag {tile_dir}

    conda 'samtools'

    input:
        val(tile_dir)

    output:
        file("tile_checked")
        
    shell:
    // Checks that read 1 and read 2 (if it exists) have the same number of reads
    '''
    if [[ $(samtools view -f 64 <(samtools cat !{tile_dir}/*.bam) -c) != $(samtools view -f 128 <(samtools cat !{tile_dir}/*.bam) -c) ]]; then
        echo "Mismatched number of reads" "!{params.flowcell} has a different number of read 1s and read 2s: !{tile_dir}"
        exit 1
    fi
    touch tile_checked
    '''
}

process determinePostRun {
    cpus 1

    input:
        val(_file_checked)


    shell:
    '''
    nextflow run !{seq_shepherd_dir}/determine_post_run.nf \
    --reads_to_transfer !{params.read_structure.count("T")} \
    --run_path !{run_path} \
    --run_uri !{params.run_uri} \
    --emails !{params.emails} \
    --flowcell !{params.flowcell} \
    --platform Illumina \
    -with-report  "!{run_path}/determinePostRun_report.html" \
    -with-timeline "!{run_path}/determinePostRun_timeline.html" \
    -with-dag "!{run_path}/determinePostRun_dag.html" \
    -with-weblog !{params.run_uri}
    '''

}

workflow {
    // We allow for multiple sample sheet naming conventions
    // Bad news if both sample sheet channels have data
    single_sample_sheets = Channel.fromPath(['./SampleSheet?csv', './*_ss.csv']) // Should be only one of these
    lane_sample_sheets = Channel.fromPath('./L*_samplesheet.csv') // Could be one per lane

    sample_sheets = single_sample_sheets.concat(lane_sample_sheets)

    // Need to get the library, barcode, etc. information
    sample_sheet_data = sample_sheets.map{it -> [it.splitCsv(header: true, skip: 23)]}.transpose().map{it -> it.flatten()}
    lanes = Channel.of(1..params.lanecount)

    // The glob in config?xml makes it not appear if the file doesn't exist
    tile_info_file = Channel.fromPath(['Data/Intensities/config?xml', 'RunInfo.xml']).first()

    writeParamsFiles(lanes.combine(sample_sheet_data))
    combineParamsFiles(writeParamsFiles.out.individual_params_files.groupTuple())
    extractIlluminaBarcodes(combineParamsFiles.out.barcode_params_file)
    findTiles(lanes.combine(tile_info_file))

    after_demux_tiles = combineParamsFiles.out.library_params_file.combine(findTiles.out.tiles.map{it -> [it[0], it[1].split("\n").toList()]}.transpose(), by: 0)
    illuminaBasecallsToSam(after_demux_tiles.combine(extractIlluminaBarcodes.out.demux_stats_file, by: 0))

    checkReadCounts(illuminaBasecallsToSam.out.tile_dir)

    if (params.demux_mode == "Internal") {
        updateReadCounts(checkReadCounts.out.collect())
        determinePostRun(checkReadCounts.out.collect())
    }
}
