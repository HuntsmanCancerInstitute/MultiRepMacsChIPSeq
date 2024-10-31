# MultiRepMacsChIPSeq - run\_normR\_enrichment

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## run\_normR\_enrichment.R

This script will run a basic enrichR function from the normR 
package to identify enriched ChIPseq regions. It requires 
an input text file with chromosome, start,and stop columns, 
along with a columns of alignment counts for both ChIP and 
reference control (Input) samples. 

Result intervals are filtered for the given Q-value FDR 
threshold, as well as a minimum read count (ChIP+Reference). 
A bedGraph file of standardized, log enrichment between 
ChIP and Reference is written for converting into a bigWig 
for visualization. Merged significant intervals for 
enrichment are written as bed files.

USAGE: 

	run_normR_enrichment.R [options]

OPTIONS:

	-i INPUT, --input=INPUT
		Input file containing count data
	
	-o OUTPUT, --output=OUTPUT
		Output file basename, default 'chip_ref' names
	
	-c CHIP, --chip=CHIP
		Name of ChIP count column
	
	-r REF, --ref=REF
		Name of reference (or Input) count column
	
	-t THRESHOLD, --threshold=THRESHOLD
		Threshold Q-value for filtering, default 0.001
	
	-m MIN, --min=MIN
		Minimum interval count sum, default 50
	
	--all
		Report all windows, not just significant
	
	-h, --help
		Show this help message and exit



