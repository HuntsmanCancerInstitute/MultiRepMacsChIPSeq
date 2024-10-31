# MultiRepMacsChIPSeq - Usage

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Examples](Examples.md)|[Applications](applications.md)|[Install](Install.md)|

## Usage guide

Below is a usage guide of the Multi-Replica Macs ChIP-Seq Wrapper. This is primarily
focused on ChIPSeq, but the same principles can (more or less) be applied for
ATAC-Seq and Cut & Run (and Cut & Tag) analyses. 

## Alignment

ChIP-Seq samples may be aligned with your favorite aligner, including NovoCraft 
[Novoalign](http://www.novocraft.com/products/novoalign/), 
[Bowtie2](http://bowtie-bio.sourceforge.net/index.shtml), or 
[BWA](https://github.com/lh3/bwa). Alignments should be converted to sorted, indexed 
Bam files using [samtools](https://github.com/samtools/samtools) or equivalent.

Duplicate alignments do not need to be removed, unless you are using Unique Molecular
Index (UMIs) tags, in which case the files should be de-duplicated first; see the
duplicate step below.

There are two metrics that are helpful to know prior to the execution of the pipeline. 
This includes estimating (single-end) or determining (paired-end) fragment sizes 
and checking duplicate levels.

### Fragment length estimation

For **single-end analysis**, fragment length must be estimated by looking at the
average offset between forward and reverse alignments flanking a peak. There are two
ways to determine this: using Macs2 `predictd` function and BioToolBox
[bam2wig](http://tjparnell.github.io/biotoolbox/apps/bam2wig.html). It's probably
recommended to run both, as each has similar but different algorithms that usually
derive similar values, and then evaluate and take the most reasonable one value. 

- Example using Macs2:

		macs2 predictd -i chip.bam -g hs --rfile chip.predictd
	
		RScript chip.predictd 

	Note that the genome size (or `hs` or `mm` shortcuts) should be provided to
	Macs2. The resulting file is an R script that will generate a PDF for
	visualization.

- Example using bam2wig:

		bam2wig.pl -i chip.bam --shift --model --out chip_shift 
	
		plot_shift_models.R --input chip_shift

	Note that bam2wig example will only sample a subset of the chromosomes unless you
	direct it to check more with `--chroms` option. It is multi-threaded with the
	`--cpu` option. It will generate two text output files, which can be processed
	into PDFs for visualization using the included
	[plot_shift_models](applications.md#plot_shift_modelsr) script as shown.

For **paired-end analysis**, only properly-mapped pairs are analyzed (singletons are 
skipped for simplicity). Knowing the mean fragment size allows you to intelligently set 
the expected peak size. Some aligners, e.g. Novoalign, will report mean fragment 
length. Otherwise, a tool such as 
[Picard CollectInsertSizeMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics)
will work.

### Duplication level determination

The pipeline is designed to explicitly handle duplicate alignments in a number of ways.
Knowing the duplication level in advance may help to set a reasonable expectation for 
how to handle duplicates. One of the methods is duplication sub-sampling; see 

You can obtain a quick determination by running the
[bam\_partial\_dedup](applications/bam_partial_dedup.md) application for each file.
If desired, capture the standard output to file, which can then be combined into a
single file with [combine\_std\_chipstats](applications/combine_std_chipstats.md). Be
sure to specify paired-end alignments with `--pe` option as necessary.

	bam_partial_dedup.pl -i file1.bam > file1.dup.txt
	
	combine_std_chipstats.pl dup_stats.txt file*.dup.txt

**Note** that it is best to exclude known sources of duplicate reads, such as 
the mitochondrial chromosome, rDNA and other high-copy number genes,repetitive 
elements, and other known sites of artificial enrichment. In the pipeline, these 
should be excluded by default, but you can manually specify the exclusion lists 
and chromosomes in `bam_partial_dedup` if available.

	bam_partial_dedup.pl -i file1.bam --chrskip chrM --exclude exclusion.bed > file1.dedup.txt

**NOTE** that if your sample has high optical duplicates, e.g. sequence from a
patterned Illumina flowcell like NovaSeq, please add the `--optical` and `--distance`
options to remove these. Optical duplicates are entirely artificial sequencing
artifacts and should not be considered. For Novaseq, set the pixel distance to at least 
2500; for NovaseqX, set the pixel distance to 200. 

## Running the pipeline

The pipeline is executed by passing all parameters to the
[multirep_macs2_pipeline](applications/multirep_macs2_pipeline.md) script. Because
of the number and complexity of the options, it's generally recommended to place the
command in a shell script and execute that. Before running, it's highly recommended
to run the command with the `--dryrun` option to check that all inputs are valid. 

The script will spawn numerous separate jobs (sub-tasks) in steps, running them in
parallel. Job control is explained below in the options section. In general, the
script is designed to run _on a single workstation or compute node_, rather than
spawning jobs through a cluster job manager. Each step command is printed to `STDOUT`
as the pipeline is run. Individual job output is captured in separate `.out.txt`
files and combined into a single log file upon completion. Progress is tracked using
a simple `output.progress.txt` file, which is used to skip finished steps when
restarting an incomplete run.

## Simple Peak call

To call peaks on a single-condition ChIP-Seq sample with multiple replicates, call
the [multirep_macs2_pipeline](applications/multirep_macs2_pipeline.md) script,
treating each replicate as a separate sample. This allows the replicates to be
compared directly. If each replicate has a separate reference, then a separate
`--control` option can be repeated for each sample in the same order.

	multirep_macs2_pipeline.pl \
	--chip file1.bam \
	--chip file2.bam \
	--chip file3.bam \
	--control input.bam \
	--name my_experiment 

Each replicate will have a separate peak call, and the replicates merged into a 
single master set of peaks for subsequent analysis and comparison.

## Multiple-sample peak calls

When multiple conditions are being tested and compared for differential binding, then
simply add additional `--chip` and `--name` arguments to the
[multirep_macs2_pipeline](applications/multirep_macs2_pipeline.md) script. In this 
case, replicates for each condition are averaged together prior to peak calling, with 
the assumption that a greater diversity and variance will be observed between the 
conditions than between the replicates within each condition. Independent replicate 
peak calls can also alternatively be made, if desired, then merged.

**NOTE** The order of multiple `--chip`, `--name`, and `--control` options is critical, 
and must be provided in the same order. If there are multiple controls, and some are 
shared between more than one ChIP but not all, then that's ok; list each control for 
each ChIP and duplicate entries will be smartly handled.

	multirep_macs2_pipeline.pl \
	--chip file1.bam,file2.bam,file3.bam \
	--chip file4.bam,file5.bam,file6.bam \
	--chip file7.bam,file8.bam \
	--control input1.bam \
	--control input2.bam \
	--control input3.bam \
	--name condition1 \ 
	--name condition2 \ 
	--name condition3 \
	--out all_chips

At the end of the pipeline, individual peaks from all conditions are merged into a 
single master list of all peaks identified. This is used in subsequent differential 
analysis to determine significant differential occupancy. A number of comparison 
plots are generated to assist in evaluation.


## Pipeline options

These are descriptions and guidance to the variety of options to the main 
[multirep_macs2_pipeline](applications/multirep_macs2_pipeline.md) script. In most 
cases, you will want to write the command in a shell script for execution due to the 
complexity. See the [examples](examples/Readme.md) folder for example scripts.

- Sample files

	There are three parameters: `--chip`, `--control`, and `-name`, for each
	sample or experimental condition. Hence, the options may be repeated, but
	the order must be maintained.
	
	The `--chip` option accepts the bam filenames. It is required. If you have
	replicates, they can be specified as a comma-delimited list (no spaces). Advanced
	users could provide a single fragment coverage bigWig file. The option is the
	same regardless	 whether your experiment is ChIP-Seq, ATAC-Seq, Cut-and-Run, or
	whatever.
	
	The `--name` option assigns the name to this sample. It is required. It's
	best to keep this name short, if possible, for legibility purposes in plots.
	
	The `--control` option accepts the bam filenames for the corresponding reference
	or Input. This parameter is optional and not used in ATAC-Seq or Cut-and-Run.
	Further, if you use the same Input for multiple ChIP experiments and is shared
	amongst all samples, it may be specified only once.
	
- Exclusion list

	Exclusion lists, also known as "black" or "gray" lists, are intervals that are
	either predetermined to be known troublesome spots ("black" list) or empirically
	determined from a reference control ("gray" list). Filtering against known hot
	spots is **highly** recommended, especially when subsampling duplicates, as these
	can be a large source of duplicate reads.
	
	By default, the pipeline will automatically generate an exclusion list from
	provided reference control (Input) Bam files, if they are present. Usually,
	empirically derived lists are superior to externally provided lists, as they
	account for duplicate and repetitive regions in the actual cells of interest. 
	
	To use a pre-determined exclusion region file, specify a BED, GFF, or
	any other file format with recognizable chromosomal coordinates.
	
		--exclude exclusion_peaks.bed \
	
	To call inherent exclusion peaks from Input Bam files, specify `input` (default).
	To disable automatic calling from exclusion peaks, specify `none`.

- Chromosome skipping

	Certain chromosomes may be skipped. Typically, alignments to the mitochondrial 
	chromosome can be a major source of duplicate and unwanted alignments that
	can skew alignment totals for genome normalization as well as total duplication
	rate. Since mitochondrial chromosomes are high in copy number (relative to 
	nuclear chromosomes) and lack chromatin structure, it's best to
	simply skip it. This is especially the case with certain techniques, such as
	ATAC-Seq, which can have extraordinarily high levels of mtDNA contamination
	because the DNA is incredibly accessible.
	
	Similar arguments can be made about other unmapped scaffolds and contigs, 
	which usually have high numbers of repetitive elements and very few genes
	(of interest). 
	
	Chromosomes can be skipped by providing a Perl regular expression. The 
	default expression, `chrM|MT|alt|Adapter|Lambda|PhiX`, includes the
	mitochondrial chromosome, alternate haplotype chromosomes, and any 
	artificial or control sequences added to a sequencing library prep (PhiX
	or lambda). However, this should be adjusted based on the particular
	genome build being used. For example, to discard unmapped scaffolds and
	contigs in human genome build hg38, the following could be used
	
		--chrskip 'chrM|chrUn|random' \
	
	The regular expression only needs to partially match, so this is fairly
	aggressive. Be careful to quote the expression in the shell script or on
	the command line to protect control characters. 
	
	**Hint**: You can quickly check your chromosome filter by simply running
	[print\_chromosome\_lengths.pl](applications/print_chromosome_lengths.md) on
	one of your bam or bigWig files - it will return sizes to standard output.

- Genome size

	The effective size of the genome in bp is required when calculating the
	background signal when determining statistical enrichment. This is the size of
	mappable space in the genome, i.e. not counting repetitive and `N` spacers.
	Pre-computed values could be used, but these are usually based on specific size
	k-mers. An alternative, empirical method is to simply determine the coverage of
	the actual provided samples, which takes into account the read length and
	sequencing depth of the given experiment. The script
	[report_mappable_space](applications/report_mappable_space.md) will do this
	automatically, unless an explicit genome size is provided.
	
	**Note** that ChIP-Seq experiments with low genomic coverage (50% or less), 
	common with early genomic sequencing technologies, will benefit more if given the 
	explicit genome size rather than empirical, as it will lead to lower expected 
	background signal and better statistical scores.

- Duplication levels

	The default behavior of the pipeline is to sub-sample duplicate alignments in each
	replicate to a consistent level, 5% (or less) by default. Alternatively, all 
	duplicates may be removed in a more traditional approach.
	
	The maximum allowed level of duplication when sub-sampling may be explicitly set.
	
		--dupfrac 0.05 \
		--optdist 2500 \
	
	Note that optical de-duplication is critical for this to work well, particularly
	when using patterned Illumina flow cells, such as Nextseq or Novaseq, which
	generates very high optical (technical) duplicates. Note that this works with
	Illumina CASAVA style read names; reads from SRA don't retain original spot
	names, so optical duplicate checking can't be performed. For an evaluation of the
	effectiveness of duplicate sub-sampling, see [De-Duplication
	Evaluation](DeDuplicationEvaluation/ReadMe.md).
	
	For a traditional approach to remove **all** duplicates, use the following settings:
	
		--dupfrac 0 \
		--maxdepth 1 \
	
	The maximum number of allowed reads at any position may be increased above 1, which
	might help with hotspots (but exclusion lists are a better approach). It is possible
	to set both sub-sampling and maximum depth.
	
	Or you may completely turn off de-duplication if you have already marked or removed
	duplicates. 
	
		--nodedup \
	
	If you are using unique molecular indexes (UMIs) or barcodes, see
	[UMIScripts](https://github.com/HuntsmanCancerInstitute/UMIScripts) as a possible 
	approach. Marked duplicate reads are **always** skipped regardless of these
	settings.
	
- Alignment filtering

	Poor quality alignments can be excluded with a minimum mapping quality, specified
	in the `MAPQ` field in the Sam file format. This is an integer, 0-255, indicating
	the confidence in uniquely mapping the alignment in the genome. Quality scores of
	0 have no confidence, single-digit are low quality, and double-digits are
	higer(er) quality. The default is to take all alignments with MAPQ of 0 or
	higher, but a minimum threshold may be set, for example
	
		--mapq 13 \
	
	As an alternative to setting a mapping quality filter, the `--fraction` option is
	available to down-weight multiple-mapping alignments, which usually have low
	mapping qualities. Rather than scoring the alignment as 1 (prior to depth
	normalization), the alignment is scored as a fraction of the number of hits,
	recorded in the alignment tag `NH` (when available), essentially as 1/`NH`. This
	allows for more comprehensive coverage while avoiding high coverage from low
	quality alignments.
	
	Secondary (Sam flag `0x100` or 256), QC fail (Sam flag `0x200` or 512), Duplicate
	(Sam flag `0x400` or 1024), and Supplementary (Sam flag `0x800` or 2048)
	alignments are always skipped. This is hard-coded and cannot be changed.
	
- Fragment size and coverage tracks

	The fragment size for single-end should be explicitly given so that all samples
	can use the same fragment size and make multiple conditions as comparable as
	possible. Use the methods described above to empirically check your fragment
	sizes from your bam files prior to running the pipeline.
	
		--size 250 \
	
	For paired-end ChIP-Seq, set the `--pe` flag and set the `--size` parameter to the
	mean observed insertion size; this value is usually reported by the alignment
	software, or you can run a utility such as 
	[Picard CollectInsertSizeMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics). 
	This value is used in generating the _lambda_ control tracks. Note that the
	fragment coverage is only derived from properly-paired alignments; singletons are
	silently discarded. If desired, the fragment size can be explicitly restricted to
	a size range with `--min` and `--max` options; sequencing depth is adjusted 
	accordingly.
	
	To expedite coverage track generation, fragment coverage tracks are binned (10 bp
	by default for ChIP). This greatly reduces computation time while minimally
	affecting coverage resolution. Lambda track binning is set automatically based on
	lambda size, but can be set manually for those obsessively inclined.
	
	All tracks are generated as Read (Fragment) Per Million (or RPM) depth-normalized 
	values. Multiple replicates are mean averaged and reported as a single track.
	
- Chromosome normalization

	For chromosome-specific normalizations, calculate the sum of alignments on the 
	chromosome of interest across all reference control (input chromatin) replicates; 
	don't use the ChIP bam as any enrichment will influence the calculation!
	Generate a normalization factor by dividing each count by the minimum count, so 
	that the condition with the lowest count will have a scaling factor of 1.0, while all
	other conditions are down-scaled to the same level. For example,
	
		--control file4.bam \
		--control file5.bam \
		--chrnorm 1,0.5678 \
		--chrapply "EBV" \

- Lambda-control

	This pipeline uses the local lambda-control background of Macs2 to account for
	chromatin bias in the region. It uses the maximum signal derived from the
	reference control track based on three size ranges: _d_ or the fragment size
	(`--size`), _small local lambda_ (`--slocal`), and _large local lambda_
	(`--llocal`). Either small or large _lambda_ may be turned off by setting the value
	to 0. Local _lambda_ can be completely turned off by using `--nolambda`. If a
	reference control bam file is not provided (generally NOT recommended, except in
	some cases such as ATAC-Seq), then a global mean is calculated from the ChIP
	fragment coverage track.

- Enrichment tracks

	Two enrichment tracks are generated: a log2 fold enrichment (log2FE) and a
	q-value (FDR) statistical enrichment track, suitable for visualization and peak
	calling. Since the fragment tracks are expressed as RPM depth-normalized, they
	are scaled automatically up to the minimum observed depth (in millions) amongst
	all of the provided bam files. This number may be overridden with the advanced
	option `--tdep` (not generally recommended).

- Peak detection

	The pipeline uses the generated q-value (FDR) track to call peaks. The
	`--cutoff` parameter is the -log10(qvalue) minimum value for calling a peak,
	so 2 is equivalent to 0.01 or 1% FDR. The minimum peak size (length) and gap
	distance for joining nearby peaks may also be explicitly set. Typically for point
	source or narrow peaks, peak size should be equivalent to the mean fragment size
	or a little larger, and gap size is around half of the fragment size, or mean
	read length. However, depending on the nature of the factor being detected and
	the efficiency of the enrichment, there can be considerable leeway in setting
	these values – do not consider these static! These three variables are probably
	the most critical in calling peaks.
	
		--cutoff 2 \
		--peaksize 250 \
		--peakgap 150 \

	When multiple samples are provided, the peaks are merged into a single Bed file, 
	representing a master list of all possible peaks identified across all conditions.
	This is used in the final analysis across all conditions for comparison purposes.

	Separate parameters are provided for broad peak calling, which are supplemental to 
	narrow peak calling, including `--broadcut` and `--broadgap`. In general, Macs2 
	uses these additional parameters to merge distant peak signals into broader 
	intervals. These are provided to the user as is, and no further analysis is 
	automatically performed with them.

- Independent replicate peak calls

	When multiple replicates are given per condition, peak calls can be made
	independently for each replicate.
	
		--independent \
	
	This is particularly helpful when the replicates do not have good correlation
	with each other or there is high background and spurious peak calls not shared
	between replicates. The replicate peaks are then merged into a single sample peak
	file. By default, the minimum number of replicates with overlapping peaks to
	accept as a final merged sample peak is `n - 1`, where `n` is the number of
	replicates. For example, if a peak is identified at gene `XYZ1` in two out of
	three replicates, it is retained, but not if it was identified in only one
	replicate. The minimum can be reset, for example to keep all peaks identified.
	This is a global value and cannot be adjusted for multiple samples.
	
		--minpeakover 1 \
	
	Mean-replicate peaks are always called in parallel with independent-replicate 
	peak calls. Analysis folders are kept separately. 

- Peak scoring

	The peaks are re-scored for their mean log2 fold enrichment and q-value score
	across all conditions. Additionally, normalized read counts are also collected
	across all replicates for use in downstream differential occupancy analysis.
	Relative spatial occupancy and enrichment are collected in bins flanking the
	peaks midpoints. The bin size can be set using `--binsize`. 
	
- Job control

	The main wrapper utilizes its own parallelization for execution, and is designed
	to run on a single work station or compute node with one or more CPU cores;
	compute cluster job management software is not needed.	Disk IO may limit
	effective throughput more so than the number of CPU cores.

	There are two options for controlling jobs and CPU usage, as some of the tools
	are multi-threaded and many of the applications can be run concurrently. The
	`--job` option indicates the number of simultaneous jobs (applications) that can
	be run concurrently. The `--cpu` option indicates the number of threads available
	to each job. The product of the two should not exceed the total number of cores
	or threads allowed for your machine. 

	When mistakes happen and the pipeline stops early (it may happen), it is usually
	possible to restart the pipeline after fixing the error. Each child-job is
	checked for pre-existing output and log files and skipped appropriately if they
	exist. 

	For advanced users, a bigWig file (fragment or lambda_control) may be specified
	for the samples instead of bam files to expedite subsequent runs.


## Evaluation

Unlike more straightforward types of analyses, ChIP-Seq is highly variable, owing to
the particular nature and behavior of the protein being analyzed. Some have nice,
narrow, tall peaks, while others are broad and diffuse. The best way to interpret
ChIP-Seq analysis is to not rely on tables of numbers, but rather to view the
generated bigWig data tracks in a genome browser, such as [IGV](https://igv.org). In
particular, load the called peak files (`.narrowPeak`), statistical enrichment
(`.qvalue.bw` bigWig files), and log2 fold enrichment (`.log2FE.bw`). Refining and
repeating the Macs2 peak calling application is not uncommon. 


## Recalling peaks with different parameters

It's not unusual, if after evaluation, to decide that the original peak calling 
parameters were not ideal and should be modified. Rather than rerunning the entire 
pipeline, peaks can now be recalled in a fraction of the time by reusing the 
existing q-value tracks using the [recall_peaks](applications/recall_peaks.md) 
script. Modifications can be made to the `--peaksize`, `--peakgap`, 
and `--cutoff` values, as well as the broad (gapped) peak parameters. Any changes 
affecting alignments, fragments, or lambda background will necessitate rerunning the 
entire pipeline.

**NOTE**: This only recalls mean-replicate peaks, and not the merged 
independent-replicate peaks generated with the `--independent` option. Recalling 
independent-replicate peaks currently requires rerunning the entire pipeline.
For example,

	recall_peaks.pl --in all_chip --out recall --dir ChIP-Seq --peaksize 150 --peakgap 100

The original pipeline `--out` value is given as the `--in` value here, and a new 
output defined. Point the directory to the same output directory as before, and the 
relevant bigWig files should be found and re-used. Peaks will be recalled, merged, and 
re-scored as before. New plots will be generated. New subfolders will be generated, 
with an incrementing digit suffix; original files will not be overwritten. The 
[intersect_peaks](applications/intersect_peaks.md) script can be used to compare 
the old and new peak calls, if desired.

## Differential peak analysis

When two or more conditions are used for ChIP, then a differential analysis can be
applied across the final merged set of peaks to identify those that significantly
differential. The R package [DESeq2](https://bioconductor.org/packages/DESeq2/) is
one recommended package, among several. The [run_DESeq2](applications/run_deseq2.md)
script is a convenient, simple, wrapper script for running such analysis. 

	run_DESeq2.R --count output_counts.txt.gz --sample output_samples.txt \
	--first chip1 --second chip2 \
	--norm --all --output differential_chip1_chip2 

The count table generated by this pipeline have already been sequencing-depth
normalized, hence the use of the `--norm` option (which sets all size factors to 1)
in the script. Since the peaks represent a tiny fraction of the genome, it is best to
use genome-wide counts for normalization, otherwise true differences between ChIP
peaks may be normalized away if DESeq2 is allowed to do it by default.

The script [generate_differential](applications/generate_differential.md) will
generate a differential track between two log2FE or coverage bigWig (or bedGraph)
files by subtracting the second from the first. Both tracks set a minimum value
before subtraction, to avoid regions that are not initially enriched. Enriched peaks
from each respective ChIP may be called by defining an absolute delta difference (no
statistical call is made).

	generate_differential.pl --in1 chip1.log2FE.bw --in2 chip2.log2FE.bw \
	--min 0.5 --delta 1 --len 200 --gap 50

Macs2 also has a differential analysis mode. This uses four input files, the fragment 
pileup and lambda control files for both ChIPs. 

	macs2 bdgdiff -t1 <file1> -t2 <file2> -c1 <file3> -c2 <file4> \
	--d1 <depth> --d2 <depth> -C <cutoff> --o-prefix <basename>

This uses log likelihood for differential detection. In my experience, the R packages 
work better.

## Manual intersection and scoring of peaks

Sometimes you want to intersect a subset of the peaks, perhaps filtering the peak
files themselves or not including a sample. To regenerate a new combined list of
peaks and re-score them for making new plots, you can follow the steps below. This is
essentially what was run automatically by the pipeline script (compare with the output 
of the pipeline). Adjust accordingly for your situation.

	intersect_peaks.pl --out new_list sample1.narrowPeak sample2.narrowPeak ...
	
	get_datasets.pl --method mean --in new_list.bed --out new_list_qvalue.txt \
	--format 3 *.qvalue.bw

	get_datasets.pl --method mean --in new_list.bed --out new_list_log2FE.txt \
	--format 3 *.log2FE.bw
	
	get_datasets.pl --method sum --in new_list.bed --out new_list_counts.txt \
	--format 0 *.count.bw
	
	cp *_samples.txt new_list_samples.txt
	
	plot_peak_figures.R -i new_list

In the above example, be sure to edit the `_samples.txt` file to reflect the new 
samples, if necessary. 

