# MultiRepMacsChIPSeq - run\_normR\_difference

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## run\_normR\_difference.R

This script will run a basic diffR function from the normR 
package to identify differential ChIPseq regions. It requires 
an input text file with chromosome, start,and stop columns, 
along with a columns of alignment counts for both ChIP1 and 
ChIP2 samples; input is not needed. 

Result intervals are filtered for the given Q-value FDR 
threshold as well as a minimum readcount (ChIP1+ChIP2). 
A bedGraph file of standardized, log difference between 
ChIP1 and ChIP2 is written for converting into a bigWig 
for visualization. Merged significant intervals for 
differential enrichment of ChIP1 (class 2) and ChIP2 
(class 1) are written as bed files.

USAGE: 

	run_normR_difference.R [options]

OPTIONS:

	-i INPUT, --input=INPUT
		Input file containing count data
	
	-o OUTPUT, --output=OUTPUT
		Output file basename, default 'first_second' names
	
	-f FIRST, --first=FIRST
		Name of first ChIP count column
	
	-s SECOND, --second=SECOND
		Name of second ChIP count column
	
	-t THRESHOLD, --threshold=THRESHOLD
		Threshold Q-value for filtering, default 0.001
	
	-m MIN, --min=MIN
		Minimum interval count sum, default 50
	
	--all
		Report all windows, not just significant
	
	-h, --help
		Show this help message and exit



