# MultiRepMacsChIPSeq - Applications

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Examples](Examples.md)|[Applications](applications.md)|[Install](Install.md)|

## Application documentation

Documentation for the applications included in this package.

- [bam\_partial\_dedup.pl](applications/bam_partial_dedup.md)

	An application for removing duplicate alignments in a bam file. Unlike
	other applications, this can sub-sample and remove a random subset
	of alignments to reach an acceptable fraction of duplication.

- [combine\_replicate\_data.pl](applications/combine_replicate_data.md)

	A simple application to combine the count data of replicates
	into single sample counts. It works with a tab-delimited text
	count table file.

- [combine\_std\_chipstats.pl](applications/combine_std_chipstats.md)

	An application to combine multiple metric values from Novocraft Novoalign
	output, [bam\_partial\_dedup](applications/bam_partial_dedup.md) duplication
	statistics, [bam2wig](http://tjparnell.github.io/biotoolbox/apps/bam2wig.html)
	shift determination, and Macs2 `predictd` shift determination.

- [generate_differential.pl](applications/generate_differential.md)

	An application to generate a differential track (bedGraph or bigWig
	format) from two track files.

- [generate\_mean\_bedGraph.pl](applications/generate_mean_bedgraph.md)

	An application to generate chromosomal mean coverage bedGraph from
	fragment coverage track file.

- [intersect_peaks.pl](applications/intersect_peaks.md)

	An application to properly intersect two or more peak files and
	generate a merged peak file along with multiple general statistics.

- [multirep\_macs2\_pipeline.pl](applications/multirep_macs2_pipeline.md)

	The main MultiRepMacsChIPSeq pipeline application.

- [peak2bed.pl](applications/peak2bed.md)

	An application to convert a `narrowPeak` file to a simple `bed` file.

- [plot\_peak\_figures.R](applications/plot_peak_figures.md)

	An R script to generate multiple plot figures from data files generated
	by the MultiRepMacsChIPSeq pipeline.

- [plot\_shift\_models.R](applications/plot_shift_models.md)

	An R script to generate a shift model plot from
	[bam2wig](http://tjparnell.github.io/biotoolbox/apps/bam2wig.html) shift
	determination files.

- [print\_chromosome\_lengths.pl](applications/print_chromosome_lengths.md)

	A simple application to generate a chromosome size file from an indexed
	track file, such as `bam` or `bigWig`.

- [recall_peaks.pl](applications/recall_peaks.md)

	An application for recalling peaks with different thresholds without having
	to re-run the entire MultiRepMacsChIPSeq pipeline.

- [report\_mappable\_space.pl](applications/report_mappable_space.md)

	An application to report the fraction of a genome that is covered by unique
	and non-unique alignments from one or more bam alignment files.

- [run_DESeq2.R](applications/run_deseq2.md)

	A simple R script to run `DESeq2` on collected peak counts, particularly
	for differential peak analysis.

- [subset_bigwig.pl](applications/subset_bigwig.md)

	A convenience utility application to subset one or more bigWig tracks to a
	single chromosome, ostensibly to download to a personal computer for evaluation
	in a genome browser.

- [update\_peak\_file.pl](applications/update_peak_file.md)

	An application to replace or update missing values in a `narrowPeak` format
	peak file using available track files.



## Pertinent Bio::ToolBox Applications

- [bam2wig.pl](http://tjparnell.github.io/biotoolbox/apps/bam2wig.html)

	An application for generating coverage tracks from bam alignment files
	in nearly every possible way imaginable.

- [data2wig.pl](http://tjparnell.github.io/biotoolbox/apps/data2wig.html)

	An application for generating `wig`, `bedGraph`, and `bigWig` files from
	a tab-delimited text file.

- [get_datasets.pl](http://tjparnell.github.io/biotoolbox/apps/get_datasets.html)

	An application for collecting data from track files over genomic features in
	a variety of methods.

- [get\_relative\_data.pl](http://tjparnell.github.io/biotoolbox/apps/get_relative_data.html)

	An application for collecting data from track files in windows flanking a
	reference point.

- [manipulate_datasets.pl](http://tjparnell.github.io/biotoolbox/apps/manipulate_datasets.html)

	An interactive command-line application for working with tab-delimited text
	data files. Faster and easier than importing into Excel to do simple row- and
	column-based manipulations.

- [manipulate_wig.pl](http://tjparnell.github.io/biotoolbox/apps/manipulate_wig.html)

	An application for performing manipulations on `wig`, `bedGraph`, and `bigWig`
	file formats.





