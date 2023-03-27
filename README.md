demux_illumina.nf
=======================

### What is this?
This is a nextflow wrapper script for [Picard's IlluminaBasecallsToSam](https://broadinstitute.github.io/picard/javadoc/picard/picard/illumina/IlluminaBasecallsToSam.html).
It is designed to work efficiently on a computing cluster, but will work on smaller machines as well.

### What does it do?
Run inside a run directory, it demultiplexes the run and produces unaligned bam (uBam) files

### To run it:
nextflow run demux_illumina.nf \
    --read_structure <read structure> \
    --flowcell <flowcell name> \
    --max_mismatches <maximum allowed mismatches> \
    --min_mismatch_delta <minimum distance to next closest barcode> \
    --max_no_calls <maximum allowed no calls (i.e., Ns)> \
    --lanecount <number of lanes> \
    --machine <sequencing machine> \
    --path_to_java <path to java> \
    --path_to_picard <path to picard>

### Configuration:
You will need to modify the path to picard variable to fit your environment.
This script is configured to use up to ~64G of RAM per lane demultiplexing job, and up to ~16G of RAM per tile bcl->uBam job.
If your machine is smaller, set the Xmx value appropriately.

### Dependencies:
- [Picard tools](https://broadinstitute.github.io/picard/) and Java runtime
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)

Patches and bug reports are welcome.
