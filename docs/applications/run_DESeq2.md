# MultiRepMacsChIPSeq - run_DESeq2

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## run_DESeq2.R

This script will run a basic DESeq2 differential analysis 
between two conditions to identify differential (or enriched 
in the case of ChIP and Input) ChIPseq regions. 

It requires an input text file with chromosome, start, and stop 
columns (or a coordinate name string), along with columns of fragment 
counts for both ChIP1 and ChIP2 (or Input reference) sample replicates.
DESeq2 requires at least 2 (ideally more) replicates per condition 
to estimate variance for normalization and significance testing.

Count data may be already normalized to a uniform genomic depth, in
which case it should indicated as such and SizeFactors will be set to 1.
Otherwise counts will be normalized by DESeq2. WARNING: Default 
normalization by DESeq2 may potentially erase any enrichment signal.

A sample condition file is required, consisting of two columns, 
unique sample identifiers and group names (conditions). If desired,
a third column could be added as an additional batch factor. Only
two groups are used in the contrast, defined with the `--first` and 
`--second` options; any additional groups and samples are ignored.

Results are filtered at the indicated threshold or alpha level 
(False Discovery Rate). Significant regions are written to a text file
and respective bed files. Log2 Fold Changes may be shrunken to reduce the
effect of low count peaks. If additional filtering on Log2FC will be 
performed, the values should be shrunk.

USAGE: 

	run_DESeq2.R [options]

OPTIONS:

	-c COUNT, --count=COUNT
		Input file containing count data

	-a SAMPLE, --sample=SAMPLE
		Sample condition file containing identifiers and conditions

	-o OUTPUT, --output=OUTPUT
		Output file basename, default 'first_second' names

	-f FIRST, --first=FIRST
		Name of first ChIP condition

	-s SECOND, --second=SECOND
		Name of second ChIP condition or Input reference

	-b, --batch
		Use third column in customized sample table as batch factor

	-t THRESHOLD, --threshold=THRESHOLD
		Threshold adjusted p-value (alpha) for significance (FDR), default 0.1

	--norm
		Input counts are already normalized

	--all
		Report all windows, not just significant

	--shrink
		Shrink Log2FC values for low count peaks

	-h, --help
		Show this help message and exit



