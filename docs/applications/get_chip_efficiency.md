# MultiRepMacsChIPSeq - get\_chip\_efficiency

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## get\_chip\_efficiency

A script to calculate the efficiency of ChIP samples. 

The signal from ChIPs is summed for all identified ChIP peaks, as well as the
remainder of the genome. The fraction of genomic signal corresponding to peaks
is reported. Typically this fraction should be several times higher for ChIPs
than for Input, and this ratio can be interpreted as the efficiency.

Alignment count (point) data is preferred, although fragment data could also be 
used (the numbers will just get very big). Note for count data that the sum of 
on- and off-peak values should be roughly equivalent to the sum of alignments 
used in the analysis, within tolerance of rounding errors for scaled datasets. 

More than one ChIP file can be provided. Each ChIP file is processed in a 
separate thread. 

VERSION: 1.2

USAGE:
	get_chip_efficiency.pl -i <peakfile> input.count.bw chip1.count.bw ...

OPTIONS:
    --in <file>             Input file of peak intervals: bed, narrowPeak, etc
    --group <file>          Optional list of groups for each chip file
    --out <file>            The output file, default input basename
    --cpu <int>             Number of threads, default 4
    --help                  Print documentation

