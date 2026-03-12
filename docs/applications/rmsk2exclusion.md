# MultiRepMacsChIPSeq - rmsk2exclusion

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## rmsk2exclusion.pl

A script to generate a new exclusion interval file based on extreme coverage
over known repetitive elements.

Repetitive elements are a known potential source of false-positive peak calls.
However, not all repetitive elements result in peaks. Most elements are
degenerate enough or have enough copies scattered across the genome that
they can be ignored by aligners and not generate significant coverage to be
a problem. But a few can generate extraordinarily high coverage that can
be called as peaks, detracting signal from real peaks. Often these are long,
complete transposable elements or centromeric elements.

While publicly available exclusion lists are available (so called "black"
lists), differences in lab strains, conditions, and aligners only make these
partially useful. Ideally, exclusions would be generated from Input chromatin
coverage, but this is not always available, particularly for ATAC-Seq,
Cut&Run, Cut&Tag, and related experiments.

This script takes a public list of known repetitive elements for a genome
(RepeatMasker files obtained from UCSC), filters for single or adjoining
elements that exceed the minimum length, scores them for normalized counts,
and identifies those with extreme coverage to be used as exclusion
intervals for bam filtering. Typically in empirical testing, the
distribution of coverage is fairly narrow, with only a tiny fraction of
elements having extreme coverage, which means having a relatively low
threshold (1 standard deviation from the mean) is sufficient for calling
outliers.

In general, all sample and replicate bam files in an experiment should be
provided. The average is taken across the samples to ensure that only those
elements with consistently high coverage are used for exclusion.

Simple differences in chromosome naming scheme (chr1 vs 1) are transparently
handled without issue; other contigs may be skipped.

VERSION: 0.1

USAGE:

	rmsk2exclusion.pl -r rmsk.txt -o exclusion.bed  File1.bam  File2.bam ...

OPTIONS:

	Input:
	  -r --rmsk <file>         Input RepeatMasker file. Required. These can often
								be obtained for your genome from
								https://genome.ucsc.edu
	  -i --input <file>        Coordinate file e.g. BED of repetitive elements to score.
								This can be used as an alternative --rmsk.
	  -o --out <file>          Output file basename. Required.
	  -d --data <file>         Specify one or more bam files to score. This option
								may be repeated for each bam file, or bam files
								may be simply appended to the end of the command.
	
	Filtering Options:
	  -l --length <int>        Minimum length of repeat element to consider (1000 bp)
	  -z --zscore <float>      Minimum z-score for scoring elements as exclusion (1)
	  --chrskip <regex>        Chromosome skip regex ($chrskip)
	  
	
	General:
	 --save                    Save the intermediate scored repeat element file.
	 -c --cpu <int>            Number of threads to use
	 --getdata <path>          Path to get_datasets.pl ()
	 --mandata <path>          Path to manipulate_datasets.pl ()
	 -h --help                 This help
	

