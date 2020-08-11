# Multi-Replica Macs ChIPSeq Wrapper

A multi-threaded wrapper for processing multi-replicate, multi-condition ChIPSeq samples

# Description

ChIPSeq and related methods, such as ATACSeq, is increasingly being used not just as 
a means of discovery ("where is my factor binding?") but also as an assay for experimental 
conditions ("how does my mutant affect factor X occupancy?"). Many traditional ChIPSeq 
programs are not necessarily well equipped to handle multiple replicates and/or sample 
conditions.

This is a wrapper application for processing ChIPSeq samples 
comprised of multiple biological replicas and/or multiple conditions in a manner to 
make comparisons as consistent and uniform as possible across samples. 

## Rationale

The venerable [Macs2](https://pypi.org/project/MACS2/) application 
provides a robust method of determining enrichment of ChIP fragments over input with a 
number of advantages: fragment-based pileup of ChIP signal versus simple counts, single 
base-pair resolution instead of sliding windows, estimation of local chromatin bias using 
multiple window sizes, and more.

However, Macs2 does not deal with replicates well (it simply adds replicates), which can 
introduce biases towards samples with greater sample depth. Comparing multiple conditions 
requires careful execution with identical parameters. Controlling duplicate reads depths to 
maximize ChIP signal, rather than discarding duplicates and ChIP signal, increases 
sensitivity between ChIP conditions.

This package aims to automate Macs2 ChIPSeq peak calling with support for multiple 
replicas and conditions while supporting newer normalization methods. Importantly, it 
will output normalized, processed bigWig enrichment files for subsequent genic analysis 
using, for example, [BioToolBox](https://github.com/tjparnell/biotoolbox) or 
[deepTools](https://deeptools.readthedocs.io/en/develop/). Several analysis and QC 
metric plots are generated after the run.


## Normalization methods

Below are some of the methods used for normalizing samples prior to peak calling.

- Same fragment length

    Running Macs2 in automatic mode with each condition can result in different 
    estimations of fragment length (assuming single-end sequencing). Since fragment 
    pileup is a basis of Macs2 calling, differences in lengths can result in differences 
    in peak intensity (longer lengths give broader and deeper peaks) and statistical 
    values. It's preferable to identify an average fragment length for all conditions and 
    use that for all conditions.

- Mean coverage of depth-normalized replicates

    Instead of simply adding replicates together, which can lead to sample bias towards 
    samples with greater sequencing depth, each sample replicate is depth-normalized 
    first and then averaged prior to generating the fragment pileups. This necessitates 
    setting a target depth for peak calling.

- Duplicate read normalization

    Removing all duplicates on the assumption that they may be PCR-derived will remove 
    biological duplication, which will especially occur at sites of strong enrichment. 
    This can destroy differences between, for example, genetic-based experimental conditions.
    Keeping all duplicates may lead to false positives due to poor library preparation 
    (true PCR duplication). The ideal solution (short of Unique Molecular Indexes or 
    barcodes) is between these two extremes: control duplication levels through 
    subsampling so that all samples have the same fraction of duplication relative to the 
    number of unique alignments. Crucial to this subsampling is first discarding optical 
    duplicates, which can artificially inflate the duplication rate. Patterned Illumina 
    flow cells, e.g. NovaSeq, has notoriously high optical duplication rates (observed 
    10-20%). This package includes an application to normalize duplicate depth. 

- Black-list interval and chromosome skipping

    Known hotspots of enrichment and repeat elements can produce false-positive peaks.
    Differences in copy-number of mitochondrial DNA or even unmapped-contigs can affect 
    depth normalization based on the total number of alignments and should be excluded. 
    This wrapper allows these regions or chromosomes to be skipped without editing your 
    bam file or altering alignment strategies.
    
- Chromosome-specific normalization

    ChIP targets on transfected vectors can have different enrichment levels compared to 
    the rest of genome or even other samples. Even targets on sex chromosomes could have 
    different enrichment levels. This wrapper allows for chromosome-specific normalization 
    factors to be applied on a per sample basis.

## Overview

Below is a general overview of the pipeline.

- Deduplicate

    If indicated, the duplicates are downsampled in all samples to the same fraction 
    using [bam_partial_dedup](applications.md#bam_partial_deduppl).

- Generate fragment coverage files

    Generate fragment coverage for the ChIP samples using
    [bam2wig](https://metacpan.org/pod/bam2wig.pl), combining depth-normalized replicates
    and applying other normalization factors, if indicated.

- Generate lambda control files

    Generate lambda control fragment coverage files using
    [bam2wig](https://metacpan.org/pod/bam2wig.pl) and Macs2 based on the maximum
    coverage of local (fragment size), small, and large lambda sizes. 

- Generate count files

    To facilitate generating count matrices later, point data (shifted start
    positions or paired fragment midpoints) count bigWig files are generated using
    [bam2wig](https://metacpan.org/pod/bam2wig.pl) for each sample replicate. By
    default, replicates are depth-normalized and scaled to the target depth (calculated 
    automatically as the minimum observed depth of all provided Bam files), or counts 
    may be left "raw" (no scaling whatsoever).

- Generate enrichment files

    Use Macs2 to generate q-value and Fold Enrichment tracks from the ChIP fragment 
    coverage and lambda control files for each ChIP condition separately.
    
- Call peaks

    Use Macs2 to call peaks from the q-value tracks for each ChIP condition separately 
    using the indicated threshold, minimum peak size, and peak gap size for merging.

- Intersect peaks

    Use BedTools to intersect the peaks from each ChIP condition into a master list
    of peaks across all ChIP conditions, as well as generate statistics of the amount
    of overlap between peaks. These can be plotted with
    [plot_peak_figures](applications.md#plot_peak_figuresr).

- Rescore peaks

	Use [get_datasets](https://metacpan.org/pod/get_datasets.pl) to generate matrices of
	log2 Fold Enrichment scores, q-value scores, and count (integer only) data for the
	master list of peaks. The log2FE and q-value scores can be plotted as heat maps with
	[plot_peak_figures](applications.md#plot_peak_figuresr). 

- Profile peaks

	Use [get_relative_data](https://metacpan.org/pod/get_relative_data.pl) to collect 
	the profile of plots using both fragment coverage and log2 Fold Enrichment scores. 
	These can be plotted as heat maps or mean line plots using 
	[plot_peak_figures](applications.md#plot_peak_figuresr). 

- Score genome

	Optionally use get_datasets to score the entire genome in windows for use in other 
	ChIPSeq analysis software, such as DESeq2 or normR (see next).
	
- Differential analysis
	
	The peaks may be evaluated for differential significance using count data and an 
	R package, such as DESeq2 or normR. This isn't part of the main pipeline per se, but 
	accessory R scripts are included which can aptly perform these functions on a per 
	sample basis. See below for more details. 

# Reference guide

Below is a rough guide on utilization.

## Alignment

ChIPSeq samples may be aligned with your favorite aligner, including NovoCraft 
[Novoalign](http://www.novocraft.com/products/novoalign/), 
[Bowtie2](http://bowtie-bio.sourceforge.net/index.shtml), or 
[BWA](https://github.com/lh3/bwa). Alignments should be converted to sorted, indexed 
bam files using [samtools](https://github.com/samtools/samtools) or equivalent.

University of Utah users using the HCI 
[pysano](https://healthcare.utah.edu/huntsmancancerinstitute/research/shared-resources/center-managed/bioinformatics/pysano/) 
system for executing jobs at [CHPC](https://www.chpc.utah.edu) can use the templates 
in the pysano folder. 

## Fragment length estimation

This may be skipped for paired-end sequenced samples.

Single-end sequenced samples should be empirically checked for fragment length. There 
are two ways to determine this: using Macs2 `predictd` function and BioToolBox 
[bam2wig](https://metacpan.org/pod/bam2wig.pl). 
It's probably recommended to run both, as each has similar but different algorithms that 
usually derive similar values, and then evaluate and take the most reasonable one value. 

- Example using Macs2:

        $ macs2 predictd -i chip.bam -g hs --rfile chip.predictd
    
        $ RScript chip.predictd 

    Note that the genome size (or `hs` or `mm` shortcuts) should be provided to Macs2. The 
    resulting file is an R script that will generate a PDF for visualization.

- Example using bam2wig:

        $ bam2wig.pl -i chip.bam --shift --model --out chip_shift 
    
        $ plot_shift_models.R --input chip_shift

	Note that bam2wig example will only sample a subset of the chromosomes unless you
	direct it to check more with `--chroms` option. It is multi-threaded with the
	`--cpu` option. It will generate two text output files, which can be processed
	into PDFs for visualization using the included
	[plot_shift_models](applications.md#plot_shift_modelsr) script as shown.

## Duplication level determination

Calculate the duplication level using the
[bam_partial_dedup](applications.md#bam_partial_deduppl) application for each
sample. If desired, capture the standard output to file, which can then be combined
into a single file with
[combine_std_chipstats](applications.md#combine_std_chipstatspl). Be sure to
specify paired-end alignments with `--pe` option as necessary.

    $ bam_partial_dedup.pl -i file1.bam > file1.dup.txt
    
    $ combine_std_chipstats.pl dup_stats.txt file*.dup.txt

**Note** that if duplicate subsampling is used (default), then it's best to avoid
sources of artificial duplicates, which can artificially inflate duplication rates.
These include the mitochondrial chromosome, rDNA and other high-copy number genes,
repetitive elements, and other known sites of artificial enrichment. If a [black
list](https://sites.google.com/site/anshulkundaje/projects/blacklists) isn't
published for your genome, you can generate one by calling peaks on the Input or
reference sample alone, and use the egregious peaks as a blacklist. You can provide
this blacklist to the [bam_partial_dedup](applications.md#bam_partial_deduppl) program.

    $ bam_partial_dedup.pl -i file1.bam --chrskip chrM --blacklist blacklist.bed > file1.dedup.txt

**NOTE** that if your sample has high optical duplicates, e.g. sequence from a
patterned Illumina flowcell like NovaSeq, please add the `--optical` and `--distance`
options to remove these. Optical duplicates are entirely artificial sequencing
artifacts and should not be considered. For NovaSeq, set the pixel distance to at least 
10,000.

	--optical --distance 10000
    
It's best to determine the level of alignment duplication in all samples, and set the
target level to the lowest observed, which is frequently the Input sample. A 5-20%
duplication rates is not uncommon in many ChIPSeq samples. If all samples have 1%
duplication or less, then de-duplication could be skipped, as it likely won't
significantly alter levels. 

Subsampling duplicate reads works best with tall, narrow ChIP enrichment, such as 
site-specific DNA-binding factors. For very broad enrichment, such as certain 
broad histone marks or certain chromatin factors, removing all duplicates may give 
better results (duplication levels may already be low).

## Simple Peak call

To call peaks on a single ChIPSeq sample with multiple replicates, call the
[multirep_macs2_pipeline](applications.md#multirep_macs2_pipelinepl) script,
giving the bam files as a comma-delimited list as shown below.

    $ multirep_macs2_pipeline.pl \
    --chip file1.bam,file2.bam,file3.bam \
    --control input.bam \
    --name my_experiment 

Using default values, this will de-duplicate the bam files to no more than 10% 
duplication level by sub-sampling, extend the single-end reads by a default 300 bp, 
average the fragment coverage depth between the replicates, generate a control-lambda track 
from the input file, generate a log2 fold-enrichment and q-value tracks, call peaks 
with a minimum threshold of 2 (1% FDR) from the q-value track, and count the number 
of reads in each peak. Many of these functions can be controlled through additional 
command-line options, for example duplication level (`--dupfrac`), fragment size 
(`--size`), and peak threshold (`--cutoff`).

The pipeline script runs a number of other applications as "jobs". The exact command 
is printed to standard output for those following at home, with the output sent to an 
output log file. At the end of the pipeline, the log files are combined into a single 
file. 

## Differential peak calls

When multiple conditions are being tested and compared for differential binding, then
simply add additional `--chip` and `--name` arguments to the
[multirep_macs2_pipeline](applications.md#multirep_macs2_pipelinepl) script. **NOTE** the
order is kept the same for each condition. If each condition has a separate reference,
then `--control` can be repeated for each as well. If there are multiple controls, and
some are shared between more than one ChIP but not all, that's ok; list each control for
each ChIP and duplicate entries will be smartly handled.

    $ multirep_macs2_pipeline.pl \
    --chip file1.bam,file2.bam,file3.bam \
    --chip file4.bam,file5.bam,file6.bam \
    --chip file7.bam,file8.bam \
    --control input.bam \
    --name condition1 \ 
    --name condition2 \ 
    --name condition3 \
    --out all_chips

At the end of the pipeline, individual peaks from all conditions are merged into a 
single master list of all peaks identified. These are then re-scored for log2 fold 
enrichment and q-value, which can be used for making heat maps for a visual 
comparison. Furthermore, read counts from each of the replicates are also collected 
for evaluating differential peaks through additional software such as 
[DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html).

### Plotting results

To prepare some clustered heat maps of the intersections and the scores across the
peaks, run the included [plot_peak_figures](applications.md#plot_peak_figuresr)
R script. It will generate a distance heat map of the jaccard intersection (spatial
overlap) between the different peak calls, as well as a heat map of the Q-value
scores and log2 Fold enrichment scores. Further, it will generate a k-means clustered
heat map of the log2 fold enrichment to identify possible interesting clusters of
differential enrichment between the conditions. Give the script the `--out` name you
provided to the wrapper:

    $ plot_peak_figures.R --input all_chips


## Pipeline options

These are descriptions and guidance to the variety of options to the main 
[multirep_macs2_pipeline](applications.md#multirep_macs2_pipelinepl) script. In most 
cases, you will want to write the command in a shell script for execution due to the 
complexity. See the Pysano folder for example scripts.

- Genome size

    The effective size of the genome in bp is required by Macs2 when calculating 
    enrichment. This is the size of mappable space in the genome, i.e. not counting 
    repetitive and `N` spacers. Default values can be provided for common species using 
    the `--species` option.

- Duplication levels

    You can set the maximum allowed level of duplication, for example to 5%.
    
        --dupfrac 0.05 \
    
    Alternatively you can set a maximum number of allowed reads at any position. This
    option might help with hotspots (but black lists are a better approach). 
    
        --dupfrac 0 \
        --maxdup 10 \
    
    For a traditional approach to remove **all** duplicates, use the following settings:
    
        --dupfrac 0 \
        --maxdup 1 \
    
    Or you may completely turn off de-duplication by using the `--nodup` flag. Use 
    this option if you have already marked duplicates, perhaps by using unique molecular 
    indexes (UMIs) or barcodes (see [UMIScripts](https://github.com/HuntsmanCancerInstitute/UMIScripts), 
    for example). Marked duplicate reads are **always** skipped regardless of these settings.
    
    In general, random subsampling is preferred over setting an arbitrary maximum depth, 
    which is analogous to cutting off all peaks at a certain height, effectively making 
    all peaks the same.
    
- Read filtering

    You can filter reads based on a minimum mapping quality, overlap with known 
    trouble hot spots (Encode blacklists or repetitive regions), or even skip entire 
    chromosomes using a Perl regular expression. Filtering against known hot spots, repetitive 
    regions, ribosomal genes, and the mitochondrial chromosome is **highly** recommended, 
    especially when subsampling duplicates, as these can be a large source of duplicate 
    reads.
    
    Secondary, supplementary, and marked duplicate reads are always skipped. 
    
        --mapq 13 \
        --blacklist hg38.blacklist.bed.gz \
        --chrskip 'chrM|MT|lambda|Adapter|PhiX' \

- Fragment size

    The fragment size for single-end should be explicitly given so that all samples 
    can use the same fragment size and make multiple conditions as comparable as possible.
    Use the methods described above to empirically check your fragment sizes from your 
    bam files prior to running the pipeline.
    
        --size 250 \
    
    For paired-end ChIPSeq, skip the size parameter and set the `--pe` flag instead. 
    This will record the span of properly-paired fragments only. Fragment size can be 
    explicitly restricted to a size range with `--min` and `--max` options.
    
    For special situations where you need to shift the coverage from the actual alignment 
    start position, use the `--shift` parameter. Provide a negative number to shift 
    "upstream" or in the 5' direction.
    
    To expedite coverage track generation, fragment coverage tracks are binned (10 bp 
    by default for ChIP). This greatly reduces computation time while minimally affecting 
    coverage resolution. Lambda track binning is set automatically based on lambda size.
    
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

    This pipeline uses the local lambda-control background of Macs2 to account for chromatin 
    bias in the region. It uses the maximum signal derived from the reference control track 
    based on three size ranges: d or the fragment size (`--size`), small local lambda 
    (`--slocal`), and large local lambda (`--llocal`). Either small or large lambda may 
    be turned off by setting the value to 0. Local lambda can be completely turned off 
    by using `--nolambda`. If a reference control bam file is not provided (NOT recommended), 
    then a chromosomal mean is calculated from the ChIP fragment coverage track.

- Enrichment tracks

    Two enrichment tracks are generated: a log2 fold enrichment (log2FE) and a q-value 
    (FDR) statistical enrichment track, suitable for visualization and peak calling. 
    Since the fragment tracks are expressed as RPM depth-normalized, they are scaled 
    automatically up to the minimum observed depth (in millions) amongst all of the 
    provided bam files. This number may be overridden with the advanced option `--tdep`.

- Peak detection

    The pipeline uses the generated q-value (FDR) track to call peaks. The `--threshold` 
    parameter is the -log10(qvalue) minimum value for calling a peak, so 2 is equivalent 
    to 0.01 or 1% FDR. The minimum peak size (length) and gap for joining nearby peaks. 
    Typically, peak size is fragment size or a little larger, and gap size is around half 
    of the fragment size.
    
        --cutoff 2 \
        --peaksize 250 \
        --peakgap 150 \
    
    Peaks are reported as simple [bed](http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
    files with 5 columns. 
    
    When multiple conditions are provided, the peaks are merged into a single Bed file, 
    representing all possible peaks identified.

- Peak scoring

    The peaks are re-scored for their mean log2 fold enrichment and q-value score across
    all conditions. Additionally, normalized read counts are also collected across all 
    replicates for use in downstream differential binding analysis.
    
- Job control

	The main wrapper utilizes its own parallelization for execution, and is designed 
	to run on a single work station or compute node with one or more CPU cores; compute 
	cluster job management software is not needed.  Disk IO may limit effective 
	throughput more so than the number of CPU cores.
	
	There are two options for controlling jobs and CPU usage, as some of the tools are 
	multi-threaded and many of the applications can be run concurrently. The `--job` 
	option indicates the number of simultaneous jobs (applications) that can be run 
	concurrently. The `--cpu` option indicates the number of threads available to each 
	job. The product of the two should not exceed the total number of cores or threads 
	allowed for your machine. 
	
	When mistakes happen and the pipeline stops early (it happens), it is usually 
	possible to restart the pipeline after fixing the error. Each child-job is checked for 
	pre-existing output and log files and skipped appropriately if they exist. 
	
	For advanced users, a bigWig file (fragment or lambda_control) may be specified for 
	the samples instead of bam files to expedite subsequent runs.

## Variation with ATAC-seq

For ATACSeq, there are two variations of analysis. 

- Fragment analysis

    An analysis of the fragments generated by ATACSeq as a measure of general DNA 
    accessibility, i.e. the higher fragment coverage indicates a higher degree of open 
    chromatin. These libraries are frequently sequenced as paired-end reads, and the 
    insert size of fragments can be restricted to specific ranges for specific analyses, 
    such as 30-120 bp for sub-nucleosomal, 121-200 bp for nucleosomal fragments, and 
    higher for multi-nucleosomal fragments. Use the `--min` and `--max` options  
    to specify acceptable paired-end size ranges. For example,
    
        --pe \
        --min 30 \
        --max 120 \
    
    For single-end, you can treat it as you would as ChIP-Seq (which would calculate 
    a mean fragment size to use), or use the cut site analysis below.
    
- Cut site analysis

    An analysis of where the cut sites are to get a higher resolution analysis of explicitly
    open DNA (or DNase Hyper Sensitive) sites. This analysis is often chosen to 
    examine potential transcription factor binding sites that frequently occur at HS sites.
    Here, paired-end information is not required and can be treated as single-end. To 
    get a pileup of signal directly over the cut site, we provide a negative shift value 
    (to shift the alignment start in the 5' direction or upstream) and then extend twice 
    the absolute shift value. For example
    
        --shift -10 \
        --size 20 \

    If de-duplicating bam files and you have paired-end alignments, use the `--deduppair`
    option to de-duplicate as paired-end alignments, but the remainder of the pipeline 
    as single-end.

# Interpretation and advanced analysis

Unlike more straightforward types of analyses, ChIPSeq is highly variable, owing to the 
particular nature and behavior of the protein being analyzed. Some have nice, narrow, 
tall peaks, while others are broad and diffuse. The best way to interpret ChIPSeq analysis 
is to not rely on tables of numbers, but rather to view the generated bigWig data tracks 
in a genome browser, such as [IGV](http://software.broadinstitute.org/software/igv/home). 
In particular, load the called peak files (`.narrowPeak`), statistical enrichment 
(`.qvalue.bw` bigWig files), and log2 fold enrichment (`.log2FE.bw`). Refining and
repeating the Macs2 peak calling application is not uncommon. 

While Macs2 generally identifies typical peaks from single, DNA-binding proteins 
reasonably well, there are instances where it breaks down. For example, very broad, weak 
peaks, such as those with histone modifications, chromatin-modifying enzymes, or even 
elongating RNA polymerase, are difficult to detect. In this case, using a window-based 
approach with another application may yield better results, as opposed to the base-pair 
level analysis that Macs2 performs. 

Differences in ChIP peaks between two or more conditions or factors is the end-goal 
for many modern experimental questions. While a k-means cluster will identify some of 
these differences in a descriptive manner, using a more rigorous approach will give a 
numerical p-value and confidence thresholds for identifying the differences.

There are two R packages recommended for these types of analysis, and both use a
count-based analysis of intervals (peaks or genomic windows),
[DESeq2](https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) and
[normR](https://bioconductor.org/packages/release/bioc/html/normr.html). The DESeq2
package uses variance between sample replicates (preferably 3) to test for
significance using a negative binomial distribution. Either differential or (ChIP vs
ChIP) or enrichment (ChIP vs Input) can be performed. The normR package uses a
binomial mixture model with expectation maximization to identify either differential
regions (ChIP vs ChIP) or enrichment (ChIP vs Input). Replicates are not used, so
replicate counts must be averaged; see the included
[combine_replicate_data](applications.md#combine_replicate_datapl) script.

## Repeating genomic search for peaks

Both DESeq2 and normR can be used. Collect counts across the genome in windows using
BioToolBox [get_datasets](https://metacpan.org/pod/get_datasets.pl) or similar.
The size of window is dependent on the nature of the peaks, but 500 bp, 1 or 2 kb may
be appropriate.

Use the [run_DESeq2](applications.md#run_deseq2r) script, specifying the Input as
the second ChIP. Use the
[run_normR_enrichment](applications.md#run_normr_enrichmentr) script, specifying
the ChIP and Input. It's best to set the cutoff for both q-value (or adjusted
p-value) as well as a minimum count; Otherwise, even low-enriched regions might get
called significant.

	run_normR_enrich.R --input chip_genome_counts.txt.gz --chip chip1 --ref input \
	--min 100 --threshold 0.001 --output chip_peaks 
	
	run_DESeq2.R --count chip_genome_counts.txt.gz --sample chip_samples.txt \
	--first chip1 --second input \
	--min 100 --threshold 0.001 --output chip_peaks 

## Differential Peak analysis

Again, both packages may be used here. Use the
[run_DESeq2](applications.md#run_deseq2r) script, specifying both ChIPs.
Alternatively, run the
[run_normR_difference](applications.md#run_normr_differencer) script. 

	run_normR_enrich.R --input chip_genome_counts.txt.gz --first chip1 --second chip2 \
	--min 100 --threshold 0.001 --output differential_chip1_chip2
	
	run_DESeq2.R --count chip_genome_counts.txt.gz --sample chip_samples.txt \
	--first chip1 --second chip2 \
	--min 100 --threshold 0.001 --output differential_chip1_chip2 

The script [generate_differential](applications.md#generate_differentialpl) will generate 
a differential track between two log2FE or coverage bigWig (or bedGraph) files by 
subtracting the second from the first. Both tracks set a minimum value before subtraction, 
to avoid regions that are not initially enriched. Enriched peaks from each respective 
ChIP may be called by defining an absolute delta difference (no statistical call is made).

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

	$ intersect_peaks.pl --out new_list sample1.narrowPeak sample2.narrowPeak ...
	
	$ get_datasets.pl --method mean --in new_list.bed --out new_list_qvalue.txt \
	--format 3 *.qvalue.bw

	$ get_datasets.pl --method mean --in new_list.bed --out new_list_log2FE.txt \
	--format 3 *.log2FE.bw
	
	$ get_datasets.pl --method sum --in new_list.bed --out new_list_counts.txt \
	--format 0 *.count.bw
	
	$ cp *_samples.txt new_list_samples.txt
	
	$ plot_figures.R -i new_list

In the above example, be sure to edit the `_samples.txt` file to reflect the new 
samples, if necessary. Also, 

## Reference genome scaling

Using a reference genome for controlling enrichment between conditions, such as 
including _Drosophila_ chromatin in your ChIP assay of human chromatin as described in 
[Orlando et al](https://doi.org/10.1016/j.celrep.2014.10.018), is an advanced technique 
of analysis. Unfortunately, this normalization scale can not be used in peak calling, 
since it artificially skews the coverage depth and breaks the assumptions of equality 
(null hypothesis) between ChIP and control required for making a statistically confident 
peak call. It is best to call peaks without normalization, preferably in a normal or wild 
type situation, and then assay those peaks with external-genome normalized coverage tracks 
in subsequent analysis. 

Earlier versions of this pipeline included options for scaling. These options remain as 
advanced options but are not detailed.

# Installation

This package includes a number of Perl and R scripts in the `bin` directory and utilizes 
several external software applications. The required external applications can be usually 
be found be searching the users' `PATH`; however, they can also be specified by using 
the appropriate command line option in the script. For convenience, it may 
be best to install everything in a [modules](https://modules.readthedocs.io/en/latest/) 
environment, a Docker image, or similar.

See the associated [Install](INSTALL.md) document for details on manually installing.

### HCI users

HCI users running the pipeline on local servers can simply load the packages into your 
environment using a module command. This can also be used with [Pysano](pysano/Readme.md).

    $ module load multirepchipseq
    $ multirep_macs2_pipeline.pl
    $ module unload multirepchipseq

### Package Installation

A standard [Module::Build](https://metacpan.org/pod/Module::Build) Perl script is 
provided for semi-automated installation of the Perl components. However, this does not 
handle external applications and dependencies. See the [Install](INSTALL.md) document 
for details.

    $ perl ./Build.PL
    $ ./Build
    $ ./Build install

# AUTHOR

    Timothy J. Parnell, PhD
    Bioinformatics Shared Resource
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

# LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. For details, see the
full text of the license in the file LICENSE.

This package is distributed in the hope that it will be useful, but it
is provided "as is" and without any express or implied warranties. For
details, see the full text of the license in the file LICENSE.
