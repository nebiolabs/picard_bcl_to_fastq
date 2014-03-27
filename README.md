ngs-picard_bcl_to_fastq
=======================

This is a wrapper script for Picard's IlluminaBasecallsToFastq. Feed it the path to a Sample names file.
it will create a new folder called "fastq" and place the demultiplexed files in there.

You'll probably need to modify the path to picard variable to fit your environment.
The defaults are for very strict bacode splitting (you may want to increase MAX_MISMATCHES and MAX_NO_CALLS).

This script works by attempting to sniff run format, barcodes and platform from the files in the run directory.
It works for GAII-x, hiseq, and miseq runs I've encountered, but YMMV.


Dependencies:
Picard tools (Java runtime)
perl
xmllint
tested in "dash" shell


Patches and bug reports are welcome.
