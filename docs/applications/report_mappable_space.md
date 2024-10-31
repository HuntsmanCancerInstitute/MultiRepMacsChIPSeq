# MultiRepMacsChIPSeq - report\_mappable\_space

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## report\_mappable\_space.pl

A script to report the mappable space of a genome based on empirical mapping 
derived from one or more bam files. For a given ChIPSeq experiment, provide all 
available Bam files.

Two numbers are reported: all mappable space reported by all alignments, 
regardless of mapping quality and number of hits to the genome, and unique 
mappability, as determined by a minimum mapping quality score (default 10).
Results are written to standard out.

USAGE:

	report_mappable_space.pl *.bam
	report_mappable_space.pl --chrskip '^chrRandom.+|chrM' *.bam

OPTIONS:

	  --in <file>         An input bam file, must be sorted and indexed
	                        Repeat for each input file
	  --qual <int>        Minimum mapping quality to be considered unique (10)
	  --chrskip <regex>   Provide a regular expression for skipping unwanted chromosomes
	  --cpu <int>         Specify the number of threads to use (4)

