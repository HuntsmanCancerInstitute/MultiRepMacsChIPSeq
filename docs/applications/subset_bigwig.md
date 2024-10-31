# MultiRepMacsChIPSeq - subset_bigwig

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## subset_bigwig.pl

A script to subset bigWig files to one specific chromosome. Useful for 
downloading a smaller file to a personal computer for visual evaluation.
New files are written in the specified directory with the same basename 
appended with the chromosome name and extension.

USAGE:

	subset_bigwig.pl -c chr1 file1.bw file2.bw ...

OPTIONS:

	-c --chrom <text>       The chromosome to subset (default chr1)
	-o --out <file>         The output directory (default ./)
	-j --job <int>          Number of simultaneous jobs, (default 2)
	--wig2bw <path>         (wigToBigWig)
	--bw2wig <path>         (bigWigToWig)
	--help                  Print documentation

