# MultiRepMacsChIPSeq - combine\_replicate\_data

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## combine\_replicate\_data.pl

A script to combine replicate values into a single value and 
written as a new data file. Replicates may be combined with a 
number of different methods. Replicate groups are specified with 
a sample file: a two-column text file with sample identifiers in 
the first column and its group identifier in the second. Sample 
identifiers must match the column name in the input replicate file. 
File compression is natively handled.

Usage: 

	combine_replicate_data.pl -i counts.txt -s samples.txt -o count_means.txt

Options:

	-i --in <file>          Input file of replicate counts
	-s --sample <file>      File of replicate samples and groups
	-m --method <text>      Method of combining: sum mean median max
							default mean
	-f --format <integer>   Number of decimals to format combined values
	-o --out <file>         Output file name


