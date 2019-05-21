# Application documentation

Documentation for the applications included in this package:

- [bam_partial_dedup.pl](#bam_partial_dedup.pl)

- [combine_replicate_data.pl](#combine_replicate_data.pl)

- [combine_std_chipstats.pl](#combine_std_chipstats.pl)

- [generate_mean_bedGraph.pl](#generate_mean_bedGraph.pl)

- [intersect_peaks.pl](#intersect_peaks.pl)

- [multirep_macs2_pipeline.pl](#multirep_macs2_pipeline.pl)

- [plot_peak_figures.R](#plot_peak_figures.R)

- [plot_shift_models.R](#plot_shift_models.R)

- [print_chromosome_lengths.pl](#print_chromosome_lengths.pl)

- [run_DESeq2.R](#run_DESeq2.R)

- [run_normR_difference.R](#run_normR_difference.R)

- [run_normR_enrichment.R](#run_normR_enrichment.R)


## bam_partial_dedup.pl

A script to remove excessive duplicate alignments down to an acceptable 
fraction. This is in contrast to traditional duplicate removers that 
remove all duplicates, retaining only one alignment per position. 

Duplicate reads may be artificial (derived from PCR duplication due to 
very low input) or natural (multiple DNA fragments enriched at few narrow 
foci, such as a ChIP peaks). Removing all duplicates can significantly 
reduce high-enrichment peaks, but not removing duplicates can lead to 
false positives. An optimal balance is therefore desirable. This is 
most important when comparing between replicates or samples.

This script has two primary modes: 
1. Randomly subsample or remove duplicate reads to reach a target 
   duplication rate. Results in a more uniform duplicate reduction 
   across the genome, which can be more typical of true PCR duplication. 
   Set the --frac and --random options below. Can also optionally set 
   the --max option to remove extreme outliers.
2. Remove duplicate reads at positions that exceed a threshold, either 
   manually set or automatically calculated, to achieve a target 
   duplication rate. Not usually recommmended as it reduces signal at  
   top peaks without addressing low level duplication across the genome.

For ChIPSeq applications, check the duplication level of the input. For 
mammalian genomes, typically 10-20% duplication is observed in sonicated 
input. ChIP samples typically have higher duplication. Set the target 
fraction of all samples to the lowest observed duplication rate.

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

Existing alignment duplicate marks (bit flag 0x400) are ignored. 

Since repetitive and high copy genomic regions are a big source of duplicate 
alignments, these regions can and should be entirely skipped by providing a 
file with recognizable coordinates. Any alignments overlapping these intervals 
are skipped in both counting and writing. 

USAGE: 
	bam_partial_dedup.pl --in in.bam
	bam_partial_dedup.pl --frac 0.xx --rand --in in.bam --out out.bam
	bam_partial_dedup.pl --m X -i in.bam -o out.bam

OPTIONS:
	--in <file>      The input bam file, should be sorted and indexed
	--out <file>     The output bam file containing unique and retained 
				     duplicates; optional if you're just checking the 
				     duplication rate.
	--pe             Bam files contain paired-end alignments and only 
				     properly paired duplicate fragments will be checked for 
				     duplication. 
	--mark           Write non-optical alignments to output and mark as 
				     duplicates with flag bit 0x400.
	--random         Randomly subsamples duplicate alignments so that the  
				     final duplication rate will match target duplication 
				     rate. Must set --frac option. 
	--frac <float>   Decimal fraction representing the target duplication 
				     rate in the final file. 
	--max <int>      Integer representing the maximum number of duplicates 
				     at each position
	--optical        Enable optical duplicate checking
	--distance       Set optical duplicate distance threshold.
				     Use 100 for unpatterned flowcell (HiSeq) or 
				     2500 for patterned flowcell (NovaSeq). Default 100.
				     Setting this value automatically sets --optical.
	--keepoptical    Keep optical duplicates in output as marked 
				     duplicates with flag bit 0x400. Optical duplicates 
				     are not differentiated from non-optical duplicates.
	--coord <string> Provide the tile:X:Y integer 1-base positions in the 
				     read name for optical checking. For Illumina CASAVA 1.8 
				     7-element names, this is 5:6:7 (default)
	--blacklist <file> Provide a bed/gff/text file of repeat regions to skip
	--chrskip <regex> Provide a regex for skipping certain chromosomes
	--seed <int>     Provide an integer to set the random seed generator to 
				     make the subsampling consistent (non-random).
	--cpu <int>      Specify the number of threads to use (4) 


## combine_replicate_data.pl

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

## combine_std_chipstats.pl

A script to combine 
- Novoalign statistics
- bam_partial_dedup (or bam_umi_dedup) statistics
- bam2wig empirical shift determination 
- Macs2 predicted shift determination

It parses this information from standard output and/or error text files from these 
programs. It will also parse the `stderr.txt` and `stdout.txt` files from Pysano job
directories. It will write out a single tab-delimited with the numbers. Sample 
names are the given input file names or directory names.

Usage: 
	Pysano directories:
		combine_std_chipstats.pl <outputfile> 1234X1/ 1234X2/ ...

	Text files:
		combine_std_chipstats.pl <outputfile> file1.txt file2.txt ...
    
## generate_mean_bedGraph.pl
  
A script to generate a chromosomal mean coverage bedGraph track 
to be used in Macs2 as the global control track when there is 
no input for generating a control_lambda chromatin bias track.
This uses a bigWig file to calculate each chromosomal mean. It 
will write out a simple bedGraph representing the genome with 
the respective mean for each chromosome.
  
Usage: `generate_mean_bedGraph.pl <file1.bw> ...`
  
It will write out a bedgraph file in the same direcotory and same 
basename appended with `.global_mean.bdg`.
  
## intersect_peaks.pl

A script to intersect two or more peak files. This is a wrapper around the bedtools 
program.

It will first merge all of the peak files into a single representative bed file.

It will then run the bedtools Jaccard statistic pairwise across all of the peak 
files and write out a merged table of the results. The Jaccard statistic measures 
the amount of spatial overlap between two peak files (intersection/union) reported 
as a fraction between 0 and 1.

Finally, it will run bedtools multiinter tool to perform a multi-way intersection 
and calculate the intervals for each category of overlap. This is parsed into a 
summary file suitable for drawing a Venn diagram based on spatial overlap.

Five files will be written:
    basename.bed                    the merged peaks
    basename.jaccard.txt            the Jaccard results in a table
    basename.n_intersection.txt     the number of intersections in a table
    basename.multi.txt              data file from multi-intersection 
    basename.spatialVenn.txt        summary of spatial overlap for each category

	USAGE: intersect_peaks.pl --out <basename> peak1.narrowPeak peak2.narrowPeak ....

	OPTIONS:
		--out basename          Provide the output basename
		--bed path              Path to bedtools (bedtools)
		--help                  Print documentation

## multirep_macs2_pipeline.pl

This is a wrapper for calling and/or comparing peaks in ChIPSeq or ATACSeq with single 
or multiple replicas using the Macs2 ChIPSeq caller. It uses BioToolBox applications to 
normalize duplicate levels and read depths between samples and replicates.

Multiple ChIP samples (experiments) may be provided by repeating the chip option as 
necessary for every experiment, factor, or antibody sample. Provide a separate name 
for each sample, in the same order.

ChIP sample replicas should be comma-delimited values to the chip option. Each 
sample could have one or more replicas. Replicas will be averaged together in a 
depth-controlled manner. If for some reason you don't want to merge replicas, then 
treat them as individual samples.

One control may be used for all samples, or sample-matched controls may be 
provided by repeating the option, keeping the same order. Control replicas may 
be provided as comma-delimited lists. If multiple, but not all, ChIP samples share 
controls, then they should still be listed individually for each ChIP; duplicate 
controls will be properly handled. If no control is available (for example, ATACSeq 
often has no genomic input), then a global mean coverage will be calculated from 
the ChIP samples and used as the control. 

Fragment size should be empirically determined by the user, especially when multiple
samples and/or replicates are being used. The same fragment size is used across all
samples and replicates to ensure equal comparisons. NOTE: even in paired-end mode, 
fragment size is used for control lambda. 

By default, this employs Macs2 local lambda chromatin-bias modeling as the reference 
track derived from the provided input. This uses three sources to model chromatin bias: 
fragment (or d in Macs2 parlance), small lambda (default 1000 bp), and 
large lambda (default 10000 bp) fragment coverage. If desired, either small or 
local lambda may be turned off by setting to 0. To completely turn off lambda, set the 
nolambda option, whereupon only the control fragment is directly used as reference. 
If no control file is provided, then the chromosomal mean from the ChIP file is used 
as a (poor) substitute. 

When using a calibration genome in the ChIPSeq (ChIP-Rx), calculate the ratio of 
target/reference alignments for each replica and sample. Provide these with the 
scale options in the same order as the bam files. Note that significant 
de-duplication levels may affect these ratios; if possible, de-duplicate first.

Advanced users may provide one processed bigWig file per ChIP or control sample. 

Options:
	 Input files
	  --chip      file1,file2...    Repeat for each sample set
	  --name      text              Repeat for each sample
	  --control   file1,file2...    Repeat if matching multiple samples
	  --chscale   number,number...  Calibration scales for each ChIP replica and set
	  --coscale   number,number...  Calibration scales for each control replica and set
 
	 Output
	  --dir       directory         Directory for writing all files (./)
	  --out       file basename     Base filename for merged output files (merged)
  
	 Genome size
	  --species   [human,mouse,fish,fly,yeast]   Default (human)
	  --genome    integer           Alternatively give effective genome size
  
	 Bam options
	  --mapq      integer           Minimum mapping quality, (0)
	  --pe                          Bam files are paired-end
	  --min       integer           Minimum paired-end size allowed (50 bp)
	  --max       integer           Maximum paired-end size allowed (500 bp)
 
	 Bam filtering options
	  --chrskip   "text"            Chromosome skip regex (chrM|MT|lambda|Adapter|PhiX)
	  --blacklist file              Bed file of repeats or hotspots to avoid
  
	 Duplication filtering
	  --nodup                       Skip deduplication
	  --dupfrac   fraction          Minimum allowed fraction of duplicates (0.1)
	  --maxdup    integer           Maximum allowed duplication depth ()
	  --optdist   integer           Maximum distance for optical duplicates (0)
									  use 100 for HiSeq, 2500 for NovaSeq
	  --savebam                     Save de-duplicated bam files

	 Fragment coverage
	  --size      integer           Predicted fragment size (single-end only, 250 bp)
	  --shift     integer           Shift the fragment, e.g. ATACSeq (0 bp)
	  --slocal    integer           Small local lambda size (1000 bp)
	  --llocal    integer           Large local lambda size (10000 bp)
	  --cbin      integer           ChIP fragment bin size (10 bp)
	  --slbin     integer           Small local lambda bin size (50 bp)
	  --llbin     integer           Large local lambda bin size (100 bp)

	 Chromosome-specific normalization
	  --chrnorm   fraction          Specific chromosome normalization factor
	  --chrapply  "text"            Apply factor to specified chromosomes
 
	 Peak calling
	  --cutoff    number            Threshold q-value for calling peaks (2) 
									 Higher numbers are more significant, -1*log10(q)
	  --tdep      integer           Average sequence depth of bam files in millions (25)
	  --peaksize  integer           Minimum peak size to call (2 x size)
	  --peakgap   integer           Maximum gap between peaks before merging (1 x size)
	  --broad                       Also perform broad (gapped) peak calling
	  --broadcut  number            Q-value cutoff for linking broad regions (1)
	  --broadgap  integer           Maximum link size between peaks in broad calls (4 x size bp)
	  --nolambda                    Skip lambda control, compare ChIP directly with control
	  --rawcounts                   Use unscaled raw counts for re-scoring peaks
	  --savebdg                     Save q-value bdg files for further custom calling
	  --window    integer           Collect counts across genome in given window size
	  --discard   number            Discard genome windows with replicate sum below number (10)
	  --repmean                     Combine replicate counts as mean for each sample set
	  --plot                        Plot figures of results
  
	 Job control
	  --cpu       integer           Number of CPUs to use per job (4)
	  --job       integer           Number of simultaneous jobs (2)

	 Application  Paths
	  --bam2wig   path             (bam2wig.pl)
	  --bamdedup  path             (bam_partial_dedup.pl)
	  --macs      path             (macs2)
	  --manwig    path             (manipulate_wig.pl)
	  --wig2bw    path             (wigToBigWig)
	  --bw2bdg    path             (bigWigToBedGraph)
	  --printchr  path             (print_chromosome_lengths.pl)
	  --data2wig  path             (data2wig.pl)
	  --getdata   path             (get_datasets.pl)
	  --meanbdg   path             (generate_mean_bedGraph.pl)
	  --bedtools  path             (bedtools)
	  --intersect path             (intersect_peaks.pl)
	  --combrep   path             (combine_replicate_data.pl)
	  --plotpeak  path             (plot_peak_figures.R)


## plot_peak_figures.R

This script generates a number of heat maps and plots for identified peak calls, 
including the following:

- heat map and cluster of jaccard (spatial overlap) statistic between peaks
- heat map and cluster of the number of peak intersctions
- heat map of the mean q-value scores for each ChIP over merged peaks
- heat map of the mean log2 fold enrichment for each ChIP over merged peaks with k-means clustering (6, 8, and 10 clusters)
- pairwise Pearson, Spearman, and Euclidean distance between all sample replicates with heat map and cluster
- PCA plot between all sample replicates

Usage: `plot_peak_figures.R [options]`
	
	Options:
		-i INPUT, --input=INPUT
			Path and basename to the multirep_macs2_pipeline combined output

		-n MIN, --min=MIN
			Minimum log2 Fold Enrichment value to plot, default -4

		-x MAX, --max=MAX
			Maximum log2 Fold Enrichment value to plot, default 4

		-q QMAX, --qmax=QMAX
			Maximum q-value to plot, default 30

		-f FORMAT, --format=FORMAT
			Format of output file: png, pdf, default png

		-h, --help
			Show this help message and exit


## plot_shift_models.R

This script will plot the predicted shift models from the BioToolBox bam2wig
application. Two plots are prepared: the mean stranded coverage around peaks, 
and the mean correlation for different shift values.

Usage: `plot_shift_models.R [options]`

Options:
	-i INPUT, --input=INPUT
		Path and basename to the bam2wig *_model.txt and *_correlations.txt files
	
	-h, --help
		Show this help message and exit


## print_chromosome_lengths.pl

A script to print out the lengths of chromosomes present in a database.
A database can be Bio::DB::SeqFeature::Store database, Bam file (`.bam`),
bigWig (`.bw`) file, bigBed (`.bb`) file, fasta (`.fa` or `.fasta`) file, or 
a directory of individual fasta files. 

Usage: `print_chromosome_lengths.pl <database>`

Options:
	-d --db "file"               Indexed database file
	-K --chrskip "text"          Chromosome skip regex
	-o --out "file"              Optional file name

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

Usage: `run_DESeq2.R [options]`

Options:
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

	-h, --help
		Show this help message and exit


## run_normR_difference.R

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

Usage: `run_normR_difference.R [options]`

Options:
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

	-h, --help
		Show this help message and exit


## run_normR_enrichment.R

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

Usage: `run_normR_enrichment.R [options]`

Options:
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

	-h, --help
		Show this help message and exit


