# Application documentation

Documentation for the applications included in this package:

- [bam_partial_dedup.pl](#bam_partial_deduppl)

- [combine_replicate_data.pl](#combine_replicate_datapl)

- [combine_std_chipstats.pl](#combine_std_chipstatspl)

- [generate_differential.pl](#generate_differentialpl)

- [generate_mean_bedGraph.pl](#generate_mean_bedgraphpl)

- [intersect_peaks.pl](#intersect_peakspl)

- [multirep_macs2_pipeline.pl](#multirep_macs2_pipelinepl)

- [plot_peak_figures.R](#plot_peak_figuresr)

- [plot_shift_models.R](#plot_shift_modelsr)

- [print_chromosome_lengths.pl](#print_chromosome_lengthspl)

- [recall_peaks.pl](#recall_peakspl)

- [report_mappable_space.pl](#report_mappable_spacepl)

- [run_DESeq2.R](#run_deseq2r)

- [run_normR_difference.R](#run_normr_differencer)

- [run_normR_enrichment.R](#run_normr_enrichmentr)

- [subset_bigwig.pl](#subset_bigwigpl)


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

This script can randomly subsample or remove duplicate reads to reach a 
target duplication rate. This results in a more uniform duplicate reduction 
across the genome, which can be more typical of true PCR duplication. 
Set the target duplication fraction rate with the --frac option below. 

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

Existing alignment duplicate marks (bit flag 0x400) are ignored. 

Since repetitive and high copy genomic regions are a big source of duplicate 
alignments, these regions can and should be entirely skipped by providing a 
file with recognizable coordinates. Any alignments overlapping these intervals 
are skipped in both counting and writing. 

Usage:
 
	bam_partial_dedup.pl --in input.bam
	bam_partial_dedup.pl --frac 0.xx --in input.bam --out output.bam

Options:

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


## generate_mean_bedGraph.pl
  
A script to generate a chromosomal mean coverage bedGraph track 
to be used in Macs2 as the global control track when there is 
no input for generating a control_lambda chromatin bias track.
This uses a bigWig or bedGraph coverage file to calculate a global 
mean. Intervals without coverage are not included in the calculation.
It will write out a simple bedGraph representing the genome 
with the respective mean for each chromosome. 

It will write out a bedgraph file in the same directory with the 
basename appended with '.global_mean.bdg'. More than one file may 
be provided at a time.

Version: 2
	
	Usage: $0 <file1.bw> ...


## intersect_peaks.pl

A script to intersect two or more peak files. This is a wrapper around the 
[BedTools](https://bedtools.readthedocs.io) program.

It will first calculate a number of descriptive statistics for the interval lengths 
for each input file, including count, sum, standard deviation, and quartiles.

It will merge all of the peak files into a single representative bed file.
Peak intervals will be renamed to the given name. 

It will run the _bedtools Jaccard_ statistic pairwise across all of the peak 
files and write out a merged table of the results. The Jaccard statistic measures 
the amount of spatial overlap between two peak files (intersection/union) reported 
as a fraction between 0 and 1.

Finally, it will run _bedtools multiinter_ tool to perform a multi-way intersection 
and calculate the intervals for each category of overlap. This is parsed into a 
summary file suitable for drawing a Venn diagram based on spatial overlap.

Six files will be written:
    basename.bed                    the merged peaks
    basename.jaccard.txt            the Jaccard results in a table
    basename.n_intersection.txt     the number of intersections in a table
    basename.multi.txt              data file from multi-intersection 
    basename.spatialVenn.txt        summary of spatial overlap for each category
    basename.lengthStats.txt        interval length statistics for each file

Usage:

	intersect_peaks.pl --out <basename> peak1.narrowPeak peak2.narrowPeak ....

Options:

		--out basename          Provide the output basename
		--name text             Provide text to rename the merged peaks
		--genome path           Provide a genome file for sort consistency
		--bed path              Path to bedtools (bedtools)
		--help                  Print documentation

## multirep_macs2_pipeline.pl

This is a wrapper for calling and/or comparing peaks in ChIPSeq or ATACSeq with
single or multiple replicas using the [Macs2](https://pypi.org/project/MACS2/)
ChIPSeq caller. It uses [BioToolBox](https://metacpan.org/pod/Bio::ToolBox)
applications to normalize duplicate levels and read depths between samples and
replicates.

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

Advanced users may provide one processed bigWig file per ChIP or control sample. 

Version: 16

Options:

	 Input files
	  --chip      file1,file2...    Repeat for each sample set
	  --name      text              Repeat for each sample
	  --control   file1,file2...    Repeat if matching multiple samples
 
	 Output
	  --dir       directory         Directory for writing all files (./)
	  --out       file basename     Base filename for merged output files (merged)
  
	 Genome size
	  --genome    integer           Specify effective mappable genome size 
	                                  (default empirically determined)
  
	 Bam options
	  --mapq      integer           Minimum mapping quality, (0)
	  --pe                          Bam files are paired-end
	  --min       integer           Minimum paired-end size allowed (50 bp)
	  --max       integer           Maximum paired-end size allowed (500 bp)
	  --fraction                    Record multiple-hit alignments as fraction of hits
 
	 Bam filtering options
	  --chrskip   "text"            Chromosome skip regex (chrM|MT|lambda|Adapter|PhiX)
	  --blacklist file              Bed file of repeats or hotspots to avoid
	                                  Determined empirically from control (Input) samples
	
	 Duplication filtering
	  --nodedup                     Skip deduplication and take everything as is
	  --dupfrac   floa              Target duplication rate for subsampling (0.05)
	  --maxdepth  integer           Maximum position alignment depth ()
	                                  set to 1 to remove all duplicates
	  --optdist   integer           Maximum distance for optical duplicates (0)
									  use 100 for HiSeq, 2500 for NovaSeq
	  --deduppair                   Run deduplication as paired-end, but coverage as single-end
	                                  e.g. for ATAC-Seq cut site analysis
	
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
	  --peaksize  integer           Minimum peak size to call (2 x fragment size)
	                                  Required for paired-end alignments.
	  --peakgap   integer           Maximum gap between peaks before merging (1 x size)
	  --broad                       Also perform broad (gapped) peak calling
	  --broadcut  number            Q-value cutoff for linking broad regions (0.5)
	  --broadgap  integer           Maximum link size between peaks in broad calls (4 x size bp)
	  --nolambda                    Skip lambda control, compare ChIP directly with control
	  --independent                 Call peaks independently for each replicate and merge
	
	 Peak scoring
	  --binsize   integer           Size of bins in 25 flanking peak bins for profile (40 bp)
	  --window    integer           Collect counts across genome in given window size
	  --discard   number            Discard genome windows with replicate sum below number (10)
	  --rawcounts                   Use unscaled raw counts for re-scoring peaks
	  --repmean                     Combine replicate counts as mean for each sample set
	  --noplot                      Do not plot figures of results
	
	 Job control
	  --cpu       integer           Number of CPUs to use per job (4)
	  --job       integer           Number of simultaneous jobs (2)
	  --dryrun                      Just print the commands without execution
	  --noorganize                  Do not organize files into subfolders when finished
	  --savebam                     Save de-duplicated bam files
	  --savebdg                     Save q-value bdg files for further custom calling

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
	  --getrel    path             (get_relative_data.pl)
	  --geteff    path             (get_chip_efficiency.pl)
	  --meanbdg   path             (generate_mean_bedGraph.pl)
	  --bedtools  path             (bedtools)
	  --intersect path             (intersect_peaks.pl)
	  --peak2bed  path             (peak2bed.pl)
	  --combrep   path             (combine_replicate_data.pl)
	  --plotpeak  path             (plot_peak_figures.R)
	  --rscript   path             (Rscript)
	  --reportmap path             (report_mappable_space.pl)


## plot_peak_figures.R

This script generates a number of heat maps and plots for identified peak calls, 
including the following:

- Scatter plot of before and after duplication counts versus non-duplicate
- PCA plot between all sample replicates based on fragment counts in peak
- Multiple pairwise correlation (Pearson, Spearman, and Euclidean distance) 
  heat maps and clusters of the counts between sample replicates
- Heat map and cluster of the number of peak intersections
- Heat map and cluster of jaccard (spatial overlap) statistic between peaks
- Pie chart of spatial overlap fraction of total merged peak coverage for 
  each sample
- Bar chart of the fraction of fragment counts in corresponding peaks 
  for each replicate, a measure of ChIP efficiency
- Heat map of the mean q-value scores for each ChIP over merged peaks
- Heat map of the mean log2 fold enrichment for each ChIP over merged peaks
  with k-means clustering (6, 8, and 10 clusters)
- Profile heat map of the fragment density over the midpoint of merged
  peaks for all samples, with and without k-means (4) clustering
- Profile heat map of the log2 Fold Enrichment over the midpoint of merged
  peaks for all samples, with and without k-means (4) clustering
- Mean profile line plot of fragment density over merged peaks
- Mean profile line plot of log2 Fold Enrichment over merged peaks

Samples are identified by color palette: Try Set1, Set2, Set3, Spectral, Dark2, 
or any other named palette in RColorBrewer. Note that excessive sample numbers 
may exceed the number the colors in a given palette (8 to 12). 

Usage: 

	plot_peak_figures.R [options]
	
Options:

		-i INPUT, --input=INPUT
			Path and basename to the multirep_macs2_pipeline combined output
		
		-n MIN, --min=MIN
			Minimum log2 Fold Enrichment value to plot, default -4
		
		-x MAX, --max=MAX
			Maximum log2 Fold Enrichment value to plot, default 4
		
		-q QMAX, --qmax=QMAX
			Maximum q-value to plot, default 30
		
		-r FMAX, --fmax=FMAX
			Maximum fragment value to plot, default 4
		
		-p --palette=NAME
			RColorBrewer palette for samples, default Set1
		
		-f FORMAT, --format=FORMAT
			Format of output file: png, pdf, default png
		
		-h, --help
			Show this help message and exit


## plot_shift_models.R

This script will plot the predicted shift models from the BioToolBox bam2wig
application. Two plots are prepared: the mean stranded coverage around peaks, 
and the mean correlation for different shift values.

Usage: 

	plot_shift_models.R [options]

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

Usage: 

	print_chromosome_lengths.pl <database>

Options:

	-d --db "file"               Indexed database file
	-K --chrskip "text"          Chromosome skip regex
	-o --out "file"              Optional file name, default db basename


## recall_peaks.pl

This script will take an existing a MultiRep ChIPSeq Pipeline results 
directory, and re-run the peak calling, merging, and rescoring steps 
under new parameters. This allows you to investigate alternative peak 
calling parameters. New QC and analysis plots are generated based on 
the new peaks. 

This does not require Bam files, nor does it re-process alignments or 
coverage files. It simply re-uses the pre-existing q-value tracks for 
making new peak calls with different parameters. If you want to change 
parameters for alignments or coverage, you will need to rerun the 
pipeline over again. Existing fragment coverage, log2 Fold Enrichment, 
and count bigWig files are used for re-scoring the new peaks. As such, 
it expects file outputs from the MultiRep ChIPSeq Pipeline. 

New peak and analysis files are placed into new subdirectories with a 
numeric suffix (allowing for multiple conditions to be run consecutively) 
without overwriting pre-existing files. Peaks from different runs can 
subsequently be compared using the intersect_peaks.pl script, if desired.

Options:

	 Input files
	  --in        file basename     Base filename for previous pipeline output
	  --dir       directory         Directory for writing all files (./)
	  --out       file basename     Base filename for new output files (merged)
  
	 Peak calling
	  --cutoff    number            Threshold q-value for calling peaks () 
	                                  Higher numbers are more significant, -1*log10(q)
	  --peaksize  integer           Minimum peak size to call ()
	  --peakgap   integer           Maximum gap between peaks before merging ()
	  --broad                       Also perform broad (gapped) peak calling
	  --broadcut  number            Q-value cutoff for linking broad regions ()
	  --broadgap  integer           Maximum link size between peaks in broad calls ()
  
	 Peak scoring
	  --binsize   integer           Size of bins in 25 flanking peak bins for profile (40)
	  --repmean                     Combine replicate counts as mean for each sample set
	  --noplot                      Do not plot figures of results
  
	 Job control
	  --cpu       integer           Number of CPUs to use per job (4)
	  --job       integer           Number of simultaneous jobs (2)
	  --dryrun                      Just print the commands without execution
	  --noorganize                  Do not organize files into subfolders when finished

	 Application  Paths
	  --macs      path             (macs2)
	  --manwig    path             (manipulate_wig.pl)
	  --wig2bw    path             (wigToBigWig)
	  --bw2bdg    path             (bigWigToBedGraph)
	  --printchr  path             (print_chromosome_lengths.pl)
	  --data2wig  path             (data2wig.pl)
	  --getdata   path             (get_datasets.pl)
	  --getrel    path             (get_relative_data.pl)
	  --geteff    path             (get_chip_efficiency.pl)
	  --meanbdg   path             (generate_mean_bedGraph.pl)
	  --bedtools  path             (bedtools)
	  --intersect path             (intersect_peaks.pl)
	  --peak2bed  path             (peak2bed.pl)
	  --combrep   path             (combine_replicate_data.pl)
	  --plotpeak  path             (plot_peak_figures.R)
	  --rscript   path             (Rscript)


## report_mappable_space.pl

A script to report the mappable space of a genome based on empirical mapping 
derived from one or more bam files. For a given ChIPSeq experiment, provide all 
available Bam files.

Two numbers are reported: all mappable space reported by all alignments, 
regardless of mapping quality and number of hits to the genome, and unique 
mappability, as determined by a minimum mapping quality score (default 10).
Results are written to standard out.

Usage:

	report_mappable_space.pl *.bam
	report_mappable_space.pl --chrskip '^chrRandom.+|chrM' *.bam

Options:

	  --in <file>         An input bam file, must be sorted and indexed
	                        Repeat for each input file
	  --qual <int>        Minimum mapping quality to be considered unique (10)
	  --chrskip <regex>   Provide a regular expression for skipping unwanted chromosomes
	  --cpu <int>         Specify the number of threads to use (4)


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

Usage: 

	run_DESeq2.R [options]

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

Usage: 

	run_normR_difference.R [options]

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
	
	--all
		Report all windows, not just significant
	
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

Usage: 

	run_normR_enrichment.R [options]

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
	
	--all
		Report all windows, not just significant
	
	-h, --help
		Show this help message and exit


## subset_bigwig.pl

A script to subset bigWig files to one specific chromosome. Useful for 
downloading a smaller file to a personal computer for visual evaluation.
New files are written in the specified directory with the same basename 
appended with the chromosome name and extension.

USAGE:

	subset_bigwig.pl -c chr1 file1.bw file2.bw ...

OPTIONS:

	-c --chrom <text>       The chromosome to subset (default chr1)
	-o --out <file>         The output directory (default ./)
	-j --job <int>          Number of simultaneous jobs, (default 2)
	--wig2bw <path>         (wigToBigWig)
	--bw2wig <path>         (bigWigToWig)
	--help                  Print documentation


