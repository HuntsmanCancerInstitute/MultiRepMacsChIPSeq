# MultiRepMacsChIPSeq - combine\_std\_chipstats

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## combine\_std\_chipstats.pl

A script to combine 
- Novoalign statistics
- bam_partial_dedup (or bam_umi_dedup) statistics
- bam2wig empirical shift determination 
- Macs2 predicted shift determination

It parses this information from standard output and/or error text files from these 
programs. It will also parse the `stderr.txt` and `stdout.txt` files from Pysano job
directories. It will write out a single tab-delimited with the numbers. Sample 
names are the given input file names or directory names.

USAGE: 

	Pysano directories:
		combine_std_chipstats.pl <outputfile> 1234X1/ 1234X2/ ...

	Text files:
		combine_std_chipstats.pl <outputfile> file1.txt file2.txt ...
    

