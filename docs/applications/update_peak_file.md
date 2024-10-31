# MultiRepMacsChIPSeq - update\_peak\_file

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## update\_peak\_file.pl

A script to fill in missing score values in a Peak file.

Running MACS2 `bdgpeakcall` or `bdgbroadcall` functions will write an
appropriate Peak output file, but will have missing score values,
specifically: 

- `signalValue`

	Fold Enrichment at summit or across broad region

- `pValue`

	P-Value at summit or across broad region

- `qValue`

	Q-Value at summit or across broad region 

One or more of these values may be updated by providing the appropriate 
score track file (currently only bigWig tracks are supported). 

Additionally, the file is automatically re-sorted with a sane chromosome 
sort, and names may be adjusted if desired, using the provided base text
that will be appended with an increasing number. Scores can also optionally 
be scaled to the UCSC-recommended range of 0-1000. 

For narrowPeak files, a summit bed file may be exported, if desired.

USAGE:

    update_peak_file.pl -i <peakfile> [options]

OPTIONS:

    -i --in <file>         Input narrowPeak file
    -o --out <file>        Output file, default input
    -n --name <text>       Base text to optionally rename the peaks
    -e --enrich <file>     The enrichment score bigWig file, e.g. log2FE
    -p --pvalue <file>     The P-Calue bigWig score file
    -q --qvalue <file>     The Q-Value bigWig score file
    -r --norm              Normalize the Score column to range 0-1000
    -u --summit            Export a summit Bed file from narrowPeak
    -h --help              Print documentation



