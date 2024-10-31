# MultiRepMacsChIPSeq - run_DESeq2

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## run_DESeq2.R

This script will run a basic DESeq2 differential analysis 
between two conditions to identify differential (or enriched 
in the case of ChIP and Input) ChIPseq regions. It requires 
an input text file with chromosome, start,and stop columns 
(or a coordinate name string), along with columns of alignment 
counts for both ChIP1 and ChIP2 (or Input reference) sample 
replicates. DESeq2 requires ideally 3 or more replicates per 
condition to estimate variance for normalization and significance 
testing. A sample condition file is required, consisting of two 
columns, sample identifiers and groups (conditions).

Result intervals are filtered for the given adjusted  
threshold as well as a minimum base count (mean of ChIP1 and 
ChIP2 replicate counts). Results are written with A bedGraph file of regularized, log2 
differences between ChIP1 and ChIP2 (or Reference) is written 
for converting into a bigWig for visualization. Merged 
significant intervals of enrichment are written as a bed file.

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

	-t THRESHOLD, --threshold=THRESHOLD
		Threshold adjusted p-value for filtering, default 0.01

	-m MIN, --min=MIN
		Minimum base count sum, default 50

	--norm
		Input counts are normalized, set SizeFactors to 1

	--all
		Report all windows, not just significant

	-h, --help
		Show this help message and exit



