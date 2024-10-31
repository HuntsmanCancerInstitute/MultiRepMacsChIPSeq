# MultiRepMacsChIPSeq - print\_chromosome\_lengths

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## print\_chromosome\_lengths.pl

A script to print out the lengths of chromosomes present in a database.
A database can be Bio::DB::SeqFeature::Store database, Bam file (`.bam`),
bigWig (`.bw`) file, bigBed (`.bb`) file, multi-fasta (`.fa` or `.fasta`) file, or 
a directory of individual fasta files. 

If more than one source file is given, then all are checked for consistency 
of chromosome names and order. An output file is only written if source 
files have the same chromosome list.

USAGE: 

	print_chromosome_lengths.pl <database1> ...
	print_chromosome_lengths.pl -K 'chrM|alt|contig|un' -o chrom.sizes *.bam

OPTIONS:

	-d --db "file"               Indexed database file, may repeat
	-K --chrskip "text"          Chromosome skip regex
	-o --out "file"              Output file name, default STDOUT

