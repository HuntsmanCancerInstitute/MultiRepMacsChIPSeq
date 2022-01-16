# Reference guide

Below is a rough guide on utilization of the Multi-Replica Macs ChIPSeq Wrapper.

## Alignment

ChIPSeq samples may be aligned with your favorite aligner, including NovoCraft 
[Novoalign](http://www.novocraft.com/products/novoalign/), 
[Bowtie2](http://bowtie-bio.sourceforge.net/index.shtml), or 
[BWA](https://github.com/lh3/bwa). Alignments should be converted to sorted, indexed 
bam files using [samtools](https://github.com/samtools/samtools) or equivalent.

University of Utah users using the HCI 
[pysano](https://healthcare.utah.edu/huntsmancancerinstitute/research/shared-resources/center-managed/bioinformatics/pysano/) 
system for executing jobs at [CHPC](https://www.chpc.utah.edu) can use the templates 
in this project's pysano folder. 

## Fragment length estimation

This is primarily for *single-end analysis only*. 

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

Generally it is safe to use a standard 5 or 10% duplication level in most standard
ChIPSeq experiments (ATACSeq and broad peaks tend to be exceptions). However, if you
wish to customize this level, you can obtain a quick determination by running the
[bam_partial_dedup](applications.md#bam_partial_deduppl) application for each sample.
If desired, capture the standard output to file, which can then be combined into a
single file with [combine_std_chipstats](applications.md#combine_std_chipstatspl). Be
sure to specify paired-end alignments with `--pe` option as necessary.

    $ bam_partial_dedup.pl -i file1.bam > file1.dup.txt
    
    $ combine_std_chipstats.pl dup_stats.txt file*.dup.txt

**Note** that it is best to exclude known sources of duplicate reads, such as 
the mitochondrial chromosome, rDNA and other high-copy number genes,repetitive 
elements, and other known sites of artificial enrichment. In the pipeline, these 
should be excluded by default, but you can manually specify the exclusion lists 
and chromosomes in `bam_partial_dedup` if available.

    $ bam_partial_dedup.pl -i file1.bam --chrskip chrM --blacklist blacklist.bed > file1.dedup.txt

**NOTE** that if your sample has high optical duplicates, e.g. sequence from a
patterned Illumina flowcell like NovaSeq, please add the `--optical` and `--distance`
options to remove these. Optical duplicates are entirely artificial sequencing
artifacts and should not be considered. For NovaSeq, set the pixel distance to at least 
2500.

	--optical --distance 2500
    
It's best to determine the level of alignment duplication in all samples, and set the
target level to the lowest observed, which is frequently the Input sample. A 5-20%
duplication rates is not uncommon in many ChIPSeq samples. If all samples have 1%
duplication or less, then de-duplication could be skipped, as it likely won't
significantly alter levels. 

Subsampling duplicate reads works best with tall, narrow ChIP enrichment, such as 
site-specific DNA-binding factors. For very broad enrichment, such as certain 
broad histone marks or certain chromatin factors, removing all duplicates may give 
better results (duplication levels may already be low).

## Running the pipeline

The pipeline is executed by passing all parameters to the
`multirep_macs2_pipeline.pl` script. Because of the number and complexity of the
options, it's generally recommended to place the command in a shell script and
execute that. Before running, it's highly recommended to run the command with the
`--dryrun` option to check that all inputs are valid. 

The script will spawn numerous separate jobs in steps, running them in parallel. Job
control is explained below in the options section. In general, the script is designed
to run on a single workstation or compute node, rather than through a cluster job
manager. Each step command is printed to `STDOUT` as the pipeline is run. Individual
job output is captured in separate `.out.txt` files and combined into a single log
file upon completion. Progress is tracked using a simple `output.progress.txt` file,
which is used to skip finished steps when restarting an incomplete run.

## Simple Peak call

To call peaks on a single-condition ChIPSeq sample with multiple replicates, call the
[multirep_macs2_pipeline](applications.md#multirep_macs2_pipelinepl) script,
treating each replicate as a separate sample. This allows the replicates to be compared 
directly. If each replicate has a separate reference, then a separate `--control` 
option can be repeated for each sample in the same order.

    $ multirep_macs2_pipeline.pl \
    --chip file1.bam \
    --chip file2.bam \
    --chip file3.bam \
    --control input.bam \
    --name my_experiment 

Each replicate will have a separate peak call, and the replicates merged into a 
single master set of peaks for subsequent analysis and comparison.

Otherwise for simplicity, replicates can simply be averaged together by using a single 
`--chip` argument and a comma-delimited list of bam files. In this case, very limited 
plots will be generated.

    $ multirep_macs2_pipeline.pl \
    --chip file1.bam,file2.bam,file3.bam \
    --control input.bam \
    --name my_experiment 

## Differential peak calls

When multiple conditions are being tested and compared for differential binding, then
simply add additional `--chip` and `--name` arguments to the
[multirep_macs2_pipeline](applications.md#multirep_macs2_pipelinepl) script. In this 
case, replicates for each condition are averaged together prior to peak calling, with 
the assumption that a greater diversity and variance will be observed between the 
conditions than between the replicates within each condition. Independent replicate 
peak calls can also alternatively be made, if desired, then merged.

**NOTE** The order of multiple `--chip`, `--name`, and `--control` options is critical, 
and must be provided in the same order. If there are multiple controls, and some are 
shared between more than one ChIP but not all, then that's ok; list each control for 
each ChIP and duplicate entries will be smartly handled.

    $ multirep_macs2_pipeline.pl \
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
[multirep_macs2_pipeline](applications.md#multirep_macs2_pipelinepl) script. In most 
cases, you will want to write the command in a shell script for execution due to the 
complexity. See the Pysano folder for example scripts.

- Exclusion list

    Exclusion lists, also known as black or gray lists, are intervals that are either 
    predetermined to be known troublesome spots (black list) or empirically determined 
    from a reference control (gray list). By default, the pipeline will automatically 
    generate an exclusion list from provided reference control (Input) Bam files. 
    Usually, empirically derived lists are superior to externally provided lists, as 
    they account for duplicate and repetitive regions in the actual cells of interest. 

- Genome size

    The effective size of the genome in bp is required when calculating the
    background signal when determining statistical enrichment. This is the size of
    mappable space in the genome, i.e. not counting repetitive and `N` spacers.
    Pre-computed values could be used, but these are usually based on specific size
    k-mers. An alternative, empirical method is to simply determine the coverage of the
    actual provided samples, which takes into account the read length and sequencing
    depth of the given experiment. The script
    [report_mappable_space](applications.md#report_mappable_spacepl) will do this
    automatically, unless an explicit genome size is provided.
    
    **Note** that ChIPSeq experiments with low genomic coverage (50% or less), 
    common with early genomic sequencing technologies, will benefit more if given the 
    explicit genome size rather than empirical, as it will lead to lower expected 
    background signal and better statistical scores.

- Duplication levels

    You can set the maximum allowed level of duplication, for example to 5%.
    
        --dupfrac 0.05 \
    
    Alternatively you can set a maximum number of allowed reads at any position. This
    option might help with hotspots (but exclusion lists are a better approach). 
    
        --dupfrac 0 \
        --maxdepth 10 \
    
    For a traditional approach to remove **all** duplicates, use the following settings:
    
        --dupfrac 0 \
        --maxdepth 1 \
    
    Or you may completely turn off de-duplication by using the `--nodedup` flag. Use
    this option if you have already marked or removed duplicates, perhaps by using 
    unique molecular indexes (UMIs) or barcodes (see
    [UMIScripts](https://github.com/HuntsmanCancerInstitute/UMIScripts), for
    example). Marked duplicate reads are **always** skipped regardless of these
    settings.
    
    In general, random subsampling is preferred over setting an arbitrary maximum depth, 
    which is analogous to cutting off all peaks at a certain height, effectively making 
    all peaks the same. 
    
    Subsampling works well for point-source datasets, such as transcription factor 
    binding or highly-restricted histone marks, such as H3K4me3, where high levels 
    of biological duplication is expected. For broad enrichment, where biological 
    duplication is rare, duplicate subsampling is less useful, and complete 
    de-duplication may be more efficient.
    
    **NOTE** Optical de-duplication is critical when sub-sampling duplicates. If your
    sequencing was done on a patterned Illumina flow cell, such as Novaseq, with
    observably high optical duplicates, you should set an optical pixel distance as
    appropriate. Note that this works with Illumina CASAVA style read names; reads
    from SRA don't retain original spot names, so optical duplicate checking can't be
    performed.
    
        --optdist 2500 \
    
- Read filtering

    You can filter reads based on a minimum mapping quality, overlap with known
    trouble hot spots (Encode blacklists or repetitive regions), or even skip entire
    chromosomes using a Perl regular expression. Filtering against known hot spots,
    repetitive regions, ribosomal genes, and the mitochondrial chromosome is
    **highly** recommended, especially when subsampling duplicates, as these can be a
    large source of duplicate reads.
    
    Secondary, supplementary, and marked duplicate reads are always skipped. 
    
        --mapq 13 \
        --blacklist hg38.blacklist.bed.gz \
        --chrskip 'chrM|MT|lambda|Adapter|PhiX' \
    
    As an alternative to setting a mapping quality filter, the `--fraction` option is
    available to down-weight multiple-mapping alignments, which usually have low
    mapping qualities. Rather than scoring the alignment as 1 (prior to depth
    normalization), the alignment is scored as a fraction of the number of hits,
    recorded in the alignment tag `NH` (when available), essentially as `1/NH`. This
    allows for more comprehensive coverage while avoiding high coverage from low
    quality alignments.

- Fragment size and coverage tracks

    The fragment size for single-end should be explicitly given so that all samples 
    can use the same fragment size and make multiple conditions as comparable as possible.
    Use the methods described above to empirically check your fragment sizes from your 
    bam files prior to running the pipeline.
    
        --size 250 \
    
    For paired-end ChIPSeq, set the `--pe` flag and set the `--size` parameter to the
    mean observed insertion size; this value is usually reported by the alignment
    software, or you can run a utility such as 
    [Picard CollectInsertSizeMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics). 
    This value is used in generating the _lambda_ control tracks. Note that the
    fragment coverage is only derived from properly-paired alignments; singletons are
    silently discarded. If desired, the fragment size can be explicitly restricted to
    a size range with `--min` and `--max` options; sequencing depth is adjusted 
    accordingly.
    
    For special situations where you need to shift the coverage from the actual
    alignment start position, use the `--shift` parameter. Provide a negative number
    to shift "upstream" or in the 5' direction. See the
    [ATACSeq](#Variation-with-ATACseq) section below for an example.
    
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
    some cases such as ATACSeq), then a global mean is calculated from the ChIP
    fragment coverage track.

- Enrichment tracks

    Two enrichment tracks are generated: a log2 fold enrichment (log2FE) and a
    q-value (FDR) statistical enrichment track, suitable for visualization and peak
    calling. Since the fragment tracks are expressed as RPM depth-normalized, they
    are scaled automatically up to the minimum observed depth (in millions) amongst
    all of the provided bam files. This number may be overridden with the advanced
    option `--tdep` (not recommended).

- Peak detection

    The pipeline uses the generated q-value (FDR) track to call peaks. The
    `--cutoff` parameter is the -log10(qvalue) minimum value for calling a peak,
    so 2 is equivalent to 0.01 or 1% FDR. The minimum peak size (length) and gap
    distance for joining nearby peaks may also be explicitly set. Typically for point
    source or narrow peaks, peak size should be equivalent to the mean fragment size
    or a little larger, and gap size is around half of the fragment size, or mean
    read length. 
    
        --cutoff 2 \
        --peaksize 250 \
        --peakgap 150 \
    
    Peaks and peak summits are reported as simple
    [bed](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) files with 5 columns. 
    
    When multiple conditions are provided, the peaks are merged into a single Bed file, 
    representing all possible peaks identified.

    When multiple replicates are given per condition, peak calls can be made
    independently for each replicate, if so desired, by providing the `--independent`
    flag. A good case for this is when the replicates do not have good correlation
    with each other and/or enrichment varies between the replicates. In this case,
    peak calling can be improved as enrichment is not averaged out between strong and
    weak replicates. The replicate peaks are merged into a single condition peak file
    prior to merging with the other conditions. 

- Peak scoring

    The peaks are re-scored for their mean log2 fold enrichment and q-value score
    across all conditions. Additionally, normalized read counts are also collected
    across all replicates for use in downstream differential occupancy analysis.
    Relative spatial occupancy and enrichment are collected in bins flanking the
    peaks midpoints. The bin size can be set using `--binsize`. 
    
    If genomic bins are desired, for example for separate peak calling or
    correlations, this option can be set with `--window` and `--discard`. 

- Job control

    The main wrapper utilizes its own parallelization for execution, and is designed
    to run on a single work station or compute node with one or more CPU cores;
    compute cluster job management software is not needed.  Disk IO may limit
    effective throughput more so than the number of CPU cores.

    There are two options for controlling jobs and CPU usage, as some of the tools
    are multi-threaded and many of the applications can be run concurrently. The
    `--job` option indicates the number of simultaneous jobs (applications) that can
    be run concurrently. The `--cpu` option indicates the number of threads available
    to each job. The product of the two should not exceed the total number of cores
    or threads allowed for your machine. 

    When mistakes happen and the pipeline stops early (it happens), it is usually
    possible to restart the pipeline after fixing the error. Each child-job is
    checked for pre-existing output and log files and skipped appropriately if they
    exist. 

    For advanced users, a bigWig file (fragment or lambda_control) may be specified
    for the samples instead of bam files to expedite subsequent runs.

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

## Re-call peaks with different parameters

It's not unusual, if after evaluation, to decide that the original peak calling 
parameters were not ideal and should be modified. Rather than re-running the entire 
pipeline, peaks can now be recalled in a fraction of the time by re-using the 
existing q-value tracks using the [recall_peaks](applications.md#recall_peakspl) 
script. Modifications can be made to the `--peaksize`, `--peakgap`, 
and `--cutoff` values, as well as the broad (gapped) peak parameters. Any changes 
affecting alignments, fragments, or lambda background will necessitate re-running the 
entire pipeline.

	recall_peaks.pl --in all_chip --out recall --dir ChIPSeq --peaksize 150 --peakgap 100

The original pipeline `--out` value is given as the `--in` value here, and a new 
output defined. Point the directory to the same output directory as before, and the 
relevant bigWig files should be found and re-used. Peaks will be recalled, merged, and 
re-scored as before. New plots will be generated. New subfolders will be generated, 
with an incrementing digit suffix; original files will not be overwritten. The 
[intersect_peaks](applications.md#intersect_peakspl) script can be used to compare 
the old and new peak calls, if desired.

## Differential peak analysis

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
samples, if necessary. 

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

