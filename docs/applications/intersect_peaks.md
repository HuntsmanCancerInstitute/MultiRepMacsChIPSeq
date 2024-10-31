# MultiRepMacsChIPSeq - intersect_peaks

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## intersect_peaks.pl

A script to intersect two or more peak interval files and generate a single 
union interval file, as well as numerous statistics describing the peaks and 
intersections. It uses the [BedTools](https://bedtools.readthedocs.io) 
application for intersection calculations.

It will calculate a number of descriptive statistics for the input files, 
including count, sum, standard deviation, and quartiles of peak lengths.

It will run a multi-way intersection with the bedtools application and parse
the output into a merged interval bed file and an overlap summary statistics
file. A minimum number of overlapping peaks may be explicitly set to restrict
the merged peaks to those represented by multiple input peak files (default 1).

It will report the pairwise number of intersections and spatial overlap 
(Jaccard statistic) between all pairwise combinations of peak intervals. 

Results may be plotted using [plot_peak_figures.R](applications/plot_peak_figures.md) 
using the out basename as input.

Seven files will be written:

    output.bed                    the merged peaks in bed format
    output.matrix.txt             boolean intersection matrix for each merged peak
    output.jaccard.txt            pairwise Jaccard statistic (bp overlap) table
    output.n_intersection.txt     pairwise intersection count table
    output.multi.txt.gz           data file from bedtools multi-intersect tool 
    output.intersection.txt       intersection statistics for each peak combination
    output.lengthStats.txt        interval length statistics for each peak input

VERSION: 6.0

USAGE:

	intersect_peaks.pl --out <basename> [options] <peak1> <peak2> ...

OPTIONS:

	  -o --out <basename>      Provide the output basename
	  -n --name <text>         Merged peak name (default output basename)
	  -m --min <integer>       Minimum number of peak overlaps to merge (default 1)
	  -a --gap <integer>       Maximum gap before merging neighboring peaks (1 bp)
	  -g --genome <path>       Provide a genome file for sort consistency
	  -b --bed <path>          Path to bedtools (bedtools)
	  -h --help                Print documentation



