# MultiRepMacsChIPSeq - annotate\_to\_tss

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## annotate\_to\_tss.pl

A script to annotate genomic intervals, such as enrichment peaks, to
neighboring gene Transcription Start Sites (TSS). 

This wraps around the `window` function of `bedtools` and collates the
output into several lists of overlapping, adjacent, and neighborhood
genes. This is predicated on the observation that regulatory regions
frequently interact with not just the closest gene but rather multiple
neighboring genes, some at considerable distances.

Note that there will be one-to-one, one-to-many, many-to-one, and many-to-many
relationships between peaks and genes. Expect duplications.

The maximum distance may be defined by the user; all TSS within this radius
will be reported. For input regions expected to be promoter-proximal, such as
H3K4me3, set this to a lower value.

This uses a transcript annotation file to extract the TSS. For best
results, a custom annotation file based on empirical expression data
should be used, preferably filtered for positive expression in the
tested samples. Otherwise, a genome annotation file will suffice.

By default, there are three output files written: one based on the input
peak files, one based on all found intersecting genes, and a relative coverage
data file.

    out.annotation.tsv          Peak table with overlapping, left, right, and
                                  neighborhood gene annotation if found
    out.genes.tsv               Gene table with overlapping, closest, adjacent,
                                  and neighborhood genes and corresponding peaks
    out.tss_profile.txt         Data table of peak spatial coverage relative to
                                  closest gene TSS

If requested, additional separate gene lists may be written:

    out.overlapping_genes.tsv   Table of all overlapping genes with peak names
    out.closest_genes.tsv       Table of the closest genes with peak names
    out.adjacent_genes.tsv      Table of all overlapping, immediate left, and 
                                  immediate right genes with peak names
    out.neighbor_genes.tsv      Table of all genes within search radius

If a gene annotation file was provided, the TSS file is also saved for future
reference.

VERSION: 0.7

USAGE:

	annotate_to_tss.pl  -i peak.bed -t GRCh38_TSS.txt

OPTIONS:

	Input:
	  -i --in <file>           Input peak file (bed, narrowPeak). Required.
	  -o --out <file>          Output file basename (default input basename)
	
	TSS Annotation (pick one only, required):
	  -t --tss <file>          TSS data file from get_gene_regions.pl
	  -a --annotation <file>   Gene annotation file (gtf, gff3, UCSC genePred)
	
	Options:
	  -d --distance <int>      Neighborhood distance for reporting in Kb (50 Kb)
	  --overlap <int>          Gap distance to consider overlapping in bp (100 bp)
	  -w --writelists          Write separate gene table lists
	
	General:
	  --bedtools <path>        Path to bedtools ()
	  --getgene <path>         Path to get_gene_regions ()
	  -h --help                Print documentation



