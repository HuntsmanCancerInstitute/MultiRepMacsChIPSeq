# MultiRepMacsChIPSeq - bam\_partial\_dedup

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## bam\_partial\_dedup.pl

A script to remove excessive duplicate alignments down to an acceptable 
fraction. This is in contrast to traditional duplicate removers that 
remove all duplicates, retaining only one alignment per position. 

Duplicate reads may be artificial (derived from PCR duplication due to 
very low input) or natural (multiple DNA fragments enriched at few narrow 
foci, such as a ChIP peaks). Removing all duplicates can significantly 
reduce high-enrichment peaks, but not removing duplicates can lead to 
false positives. An optimal balance is therefore desirable. This is 
most important when comparing between replicates or samples.

This script can randomly subsample or remove duplicate reads to reach a 
target duplication rate. This results in a more uniform duplicate reduction 
across the genome, which can be more typical of true PCR duplication. 
Set the target duplication fraction rate with the `--frac` option below. 

This script can also simply remove excessive duplicate reads at positions 
that exceed a specified target threshold. This can be set either alone or
in combination with the random subsample. Using alone is generally not 
recommmended, as it reduces signal at extreme peaks without addressing 
low level duplication elsewhere across the genome.

For ChIPSeq applications, check the duplication level of the input. For 
mammalian genomes, typically 5-20% duplication is observed in sonicated 
input. For very strong enrichment of certain targets, it's not unusual 
to see higher duplication rates in ChIP samples than in Input samples. 
Generally, set the target fraction of all samples to the lowest observed 
duplication rate.

Single-end aligment duplicates are checked for start position, strand, and 
calculated alignment end position to check for duplicates. Because of this, 
the numbers may be slightly different than calculated by traditional duplicate 
removers.

Paired-end alignments are treated as fragments. Only properly paired 
alignments are considered; singletons are skipped. Fragments 
are checked for start position and fragment length (paired insertion size) 
for duplicates. Random subsampling should not result in broken pairs.

Optical duplicates, arising from neighboring clusters on a sequencing flow 
cell with identical sequence, may now be checked. When random subsampling 
duplicates, optical duplicates should critically be ignored. This is highly 
recommended for patterned flow cells from Illumina NovaSeq or NextSeq. Set a 
distance of 100 pixels for unpatterned (Illumina HiSeq) or at least 10000 for  
patterned (NovaSeq). By default, optical duplicate alignments are not written 
to output. To ONLY filter for optical duplicates, set `--max` to a very high number.
Note that tile-edge duplicates are not counted as such.

Existing alignment duplicate marks (bit flag `0x400`) are ignored. 

Since repetitive and high copy genomic regions are a big source of duplicate 
alignments, these regions can and should be entirely skipped by providing a 
file with recognizable coordinates. Any alignments overlapping these intervals 
are skipped in both counting and writing. 

USAGE:

	bam_partial_dedup.pl --in input.bam
	bam_partial_dedup.pl --frac 0.xx --in input.bam --out output.bam

OPTIONS:

	--in <file>        The input bam file, should be sorted and indexed
	--out <file>       The output bam file containing unique and retained 
				         duplicates; optional if you're just checking the 
				         duplication rate.
	--pe               Bam files contain paired-end alignments and only 
				         properly paired duplicate fragments will be checked for 
				         duplication. Singletons are silently dropped.
	--qual <int>       Skip alignments below indicated mapping quality (0)
	--mark             Write non-optical alignments to output and mark as 
				         duplicates with flag bit 0x400.
	--frac <float>     Decimal fraction representing the target duplication
				         rate in the final file. 
	--max <int>        Integer representing the maximum number of alignments 
				         at each position. Set to 1 to remove all duplicates.
	--optical          Enable optical duplicate checking
	--distance         Set optical duplicate distance threshold.
				         Use 100 for unpatterned flowcell (HiSeq) or 
				         2500 for patterned flowcell (NovaSeq). Default 100.
				         Setting this value automatically sets --optical.
	--report           Write duplicate distance report files only, no de-duplication
	--keepoptical      Keep optical duplicates in output as marked 
				        duplicates with flag bit 0x400. Optical duplicates 
				         are not differentiated from non-optical duplicates.
	--coord <string>   Provide the tile:X:Y integer 1-base positions in the 
				         read name for optical checking. For Illumina CASAVA 1.8 
				         7-element names, this is 5:6:7 (default)
	--blacklist <file> Provide a bed/gff/text file of repeat regions to skip
	--chrskip <regex>  Provide a regex for skipping certain chromosomes
	--seed <int>       Provide an integer to set the random seed generator to 
				         make the subsampling consistent (non-random).
	--cpu <int>        Specify the number of threads to use (4) 
	--verbose          Print more information
	--help             Print full documentation



