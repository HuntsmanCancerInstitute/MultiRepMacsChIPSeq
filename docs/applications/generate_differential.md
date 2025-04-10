# MultiRepMacsChIPSeq - generate\_differential

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## generate_differential.pl

A script to generate differential enrichment track files. These are 
generated from two provided bedGraph or bigWig tracks. These may be 
fold enrichment (e.g. log2FE) or fragment coverage tracks. This works 
particularly well with nucleosome fragment coverage to identify 
changes in nucleosome density.

Input files are first "zeroed", where any negative scores are reset 
to zero (or values set to a minimum defined value), to avoid examining 
regions with little or no initial enrichment. A differential track is 
then generated by subtracting input2 from input1. Complimentary 
enrichment tracks are then generated from the differential track, 
"positive" being enriched for input1 and "negative" being enriched 
for input2. The original differential track may be kept if desired.

New peaks may be called from these respective enrichment tracks, 
if so desired. If new peaks are to be re-called, specify the minimum 
delta, peak length, and gap length. 

NOTE: Peaks are called on given absolute thresholds and NOT statistical 
confidence. Care should be taken to use a reasonably confident threshold.

USAGE: 

	generate_differential.pl -1 <chip1.log2FE.bw> -2 <chip2.log2FE.bw>

OPTIONS:

	  Required 
		-1 --in1 <file>         The bw or bdg file for ChIP-1
		-2 --in2 <file>         The bw or bdg file for ChIP-2
	
	  Differential
		--min <float>           Set the minimum value to keep in input (0)
	    --scale <float>         If necessary, optional scaling factor for RPM files
	                               only used for converting input, not output
		--keepdiff              Keep the differential file
	
	  Re-call peaks (optional)
		--delta <float>         Threshold delta score to call a peak
		--len <int>             Minimum length of peak in bp
		--gap <int>             Minimum length of gap in bp
	
	  Paths
		--macs <path>           ($macs)
		--manwig <path>         ($manwig)
		--w2bw <path>           ($wig2bw)
		--bw2w <path>           ($bw2bdg)
	
	  General
		--outdir <path>         Alternate output directory
		--bw                    Convert output differential files to bw
								   default true if input is bw
								   use --nobw for false
		--db <file>             Indexed database for converting to bw
								   default uses input bigWig
		--help                  Print documentation


