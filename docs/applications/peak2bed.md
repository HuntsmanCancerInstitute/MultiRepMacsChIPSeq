# MultiRepMacsChIPSeq - peak2bed

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## peak2bed.pl

A script to convert Macs2 narrowPeak and gappedPeak files into simpler bed files. 
Specifically, discard unused p-value and q-value columns after running Macs2 
bdgpeakcall or bdgbroadcall. 

This script will also resort the peaks properly (numerically, as opposed to 
alphabetically) and rename the peaks. Score values may be normalized to a 
standard range of 0-1000. 

For narrowPeak files, one or two files may be written:

    a 5-column BED file of peak intervals, including original score value
    a 4-column BED file of the summit position 

Files will use the same path and base name, but change the extension to 
`.bed` and `.summit.bed`, respectively.

For gappedPeak files, only one file will be written, a 12-column bed file, 
with extension `.gapped.bed`.

Version: 2

USAGE:

	peak2bed.pl peak1.narrowPeak [peak1.gappedPeak] ....

OPTIONS

    -i --input  <file>   Specify the input file
    -o --output <file>   Specify the output file (default input basename)
    -n --name   <text>   Specify the name text (default input basename)
    -s --summit          write a summit bed file (narrowPeak input only)
                            default true, use --nosummit to turn off
    -r --norm            Normalize the Score column to range 0-1000
    -h --help            print this help
    
    Multiple input and ouput files and names may be specified via options,
    simply repeat as necessary.

