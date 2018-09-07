# Multi-Replica Macs ChIPSeq Wrapper

A multi-threaded wrapper for processing multi-replicate, multi-condition ChIPSeq samples

# Description

ChIPSeq and related methods, such as ATACSeq, is increasingly being used not just as 
a means of discovery ("where is my factor binding?") but also as an assay for experimental 
conditions ("how does my mutant affect factor X occupancy?"). Many ChIPSeq programs are 
not necessarily well equipped to handle multiple replicates and/or sample conditions.

This is a wrapper application for processing ChIPSeq samples 
comprised of multiple biological replicas and/or multiple conditions in a manner to 
make comparisons as consistent and uniform as possible across samples. 

## Rationale

The venerable [Macs2](https://pypi.org/project/MACS2/2.1.1.20160309/) application 
provides a robust method of determining enrichment of ChIP fragments over input with a 
number of advantages: fragment-based pileup of ChIP signal versus simple counts, single 
base-pair resolution instead of sliding windows, estimation of local chromatin bias using 
multiple window sizes, and more.

However, Macs2 does not deal with replicates well (it simply adds replicates), which can 
introduce biases towards samples with greater sample depth. Comparing multiple conditions 
requires careful execution with identical parameters. New techniques for normalization, 
such as reference genomes, are not supported.

This package aims to automate Macs2 ChIPSeq peak calling with support for multiple 
replicas and conditions while supporting newer normalization methods.

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
    first and then averaged prior to generating the fragment pileups.

- Duplicate read normalization

    Removing all duplicates on the assumption that they may be PCR-derived will remove 
    biological duplication, which will especially occur at sites of strong enrichment. 
    This can destroy differences between, for example, genetic-based experimental conditions.
    Keeping all duplicates may lead to false positives due to poor library preparation 
    (true PCR duplication). The ideal solution (short of Unique Molecular Indexes or 
    barcodes) is between these two extremes: control duplication levels through 
    subsampling so that all samples have the same fraction of duplication relative to the 
    number of unique alignments. This package includes an application to normalize 
    duplicate depth.

- Black-list interval and chromosome skipping

    Known hotspots of enrichment and repeat elements can produce false-positive peaks.
    Differences in copy-number of mitochondrial DNA or even unmapped-contigs can affect 
    depth normalization based on the total number of alignments and should be excluded. 
    This wrapper allows these regions or chromosomes to be skipped without editing your 
    bam file or altering alignment strategies.
    
- Reference genome normalization

    Including a reference genome in your ChIP samples, for example mixing reference Drosophila 
    chromatin with human chromatin, allows one to normalize the pull-down efficiency and 
    compensate for increased or decreased affinity. This wrapper includes genome-wide 
    normalization factors on a per sample basis.

- Chromosome-specific normalization

    ChIP targets on transfected vectors can have different enrichment levels compared to 
    the rest of genome or even other samples. Even targets on sex chromosomes could have 
    different enrichment levels. This wrapper allows for chromosome-specific normalization 
    factors to be applied on a per sample basis.

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

	Note that bam2wig example will only sample a subset of the chromosomes unless you direct 
	it to check more with `--chroms` option. It is multi-threaded with the `--cpu` option. 
	It will generate two text output files, which can be processed into PDFs for 
	visualization using the included `plot_shift_models.R` script as shown.

## Duplication level determination

It's best to determine the level of alignment duplication in all samples, and set the 
target level to the lowest observed. Typically 5-20% duplication rates is not uncommon in 
many ChIPSeq samples. If all samples have 1-2% duplication or less, then de-duplication 
could be skipped, as it likely won't significantly alter levels. 

    $ bam_partial_dedup.pl -i file1.bam 

## Sample Peak call

To call peaks on a single ChIPSeq sample with multiple replicates, call the 
`multirep_macs2_pipeline.pl` script, giving the bam files as a comma-delimited list 
as shown below.

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

## Differential peak calls

When multiple conditions are being tested and compared for differential binding, then 
simply add additional `--chip` and `--name` arguments. If each condition has a separate 
reference, then `--control` can be repeated for each as well.

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
enrichment and q-value scored, which can be used for making heat maps for a visual 
comparison. Furthermore, read counts from each of the replicates are also collected 
for evaluating differential peaks through additional software such as 
[DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html).

To prepare some clustered heat maps of the intersections and the scores across the 
peaks, run the included `plot_figures.R` R script. It will generate a distance heat map 
of the jaccard intersection (spatial overlap) between the different peak calls, as well 
as a heat map of the Q-value scores and log2 Fold enrichment scores. Further, it will 
generate a k-means clustered heat map of the log2 fold enrichment to identify possible 
interesting clusters of differential enrichment between the conditions. Give the script 
the `--out` name you provided to the wrapper:

    $ plot_peak_figures.R --input all_chips


## Pipeline options

- Duplication levels

    You can set the maximum allowed level of duplication, for example to 5%.
    
        --dupfrac 0.05 \
    
    Alternatively you can set a maximum number of allowed reads at any position.
    
        --dupfrac 0 \
        --maxdup 10 \
    
    Or you may completely turn off de-duplication by using the `--nodup` flag. Use 
    this option if you have already marked duplicates, perhaps by unique molecular 
    indexes (UMIs) or barcodes.
    
- Read filtering

    You can filter reads based on a minimum mapping quality, overlap with known 
    trouble hot spots (Encode blacklists or repetitive regions), or even skip entire 
    chromosomes using a Perl regular expression. Filtering against known hot spots, repetitive 
    regions, ribosomal genes, and the mitochondrial chromosome is *highly* recommended, 
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
    This will record the span of properly-paired fragments only.
    
    For special situations where you need to shift the coverage from the actual alignment 
    start position, use the `--shift` parameter. Provide a negative number to shift 
    "upstream" or in the 3' direction.
    
    To expedite coverage track generation, fragment coverage tracks are binned (10 bp 
    by default for ChIP). This greatly reduces computation time while minimally affecting 
    coverage resolution. Lambda track binning is set automatically based on lambda size.
    
    All tracks are generated as Read Per Million (RPM) depth-normalized values. Multiple 
    replicates are mean averaged and reported as a single track.

- Genome or chromosome normalization

    When normalizing for reference genomes, such as including Drosophila chromatin in 
    your ChIP assay of human chromatin, include the normalization factor. This can be 
    calculated by summing the number of alignments to the reference genome, and calculate 
    with the following formula:
    
        factor = 1 / (reference_count / 1000000)
    
    This factor will normalize the bam coverage as reads per million reference-genome 
    mapped. It is essential to do this for _both_ ChIP and Input; shared input for 
    multiple biological replicates is not recommended. Provide the normalization factor 
    for each replicate. 
    
        --chip file1.bam,file2.bam,file3.bam \
        --chscale 0.692933,1.820191,0.814981 \
        --control file4.bam, file5.bam, file6.bam \
        --coscale 0.262140,0.301934,0.286234 \
    
    For chromosome-specific normalizations, calculate the sum of alignments on the 
    chromosome of interest across all reference control (input chromatin) replicates; 
    don't use the ChIP bam as any enrichment will influence the calculation!  
    Generate a normalization factor by dividing each count by the minimum count, so 
    that the condition with the lowest count will have a scaling factor of 1.0, while all 
    other conditions are down-scaled to the same level. For example,
    
        --control file4.bam \
        --control file5.bam \
        --chrnorm 1,0.5678 \
        --chrapply "EBV"

- Lambda-control

    This pipeline uses the local lambda-control background of Macs2 to account for chromatin 
    bias in the region. It uses the maximum signal derived from the reference control track 
    based on three size ranges: d or the fragment size (`--size`), small local lambda 
    (`--slocal`), and large local lambda (`--llocal`). Either small or local lambda may 
    be turned off by setting the value to 0. Local lambda can be completely turned off 
    by using `--nolambda`. If a reference control bam file is not provided (NOT recommended), 
    then a chromosomal mean is calculated from the ChIP fragment coverage track.

- Enrichment tracks

    Two enrichment tracks are generated: a log2 fold enrichment (log2FE) and a q-value 
    (FDR) statistical enrichment track, suitable for visualization and peak calling. 
    Since the fragment tracks are expressed as RPM depth-normalized, the number of 
    alignments from the replicates should be reported as the desired target depth. 
    Ideally, this should represent the minimum number of alignments across all replicates 
    and conditions, as down-scaling higher depths is preferable over up-scaling very 
    low sequence depths. Express this as a rounded number of millions of reads. For 
    example, if the minimum depth is 16,768,123 reads, set target depth as
    
        --tdep 17

- Peak detection

    The pipeline uses the generated q-value (FDR) track to call peaks. The `--threshold` 
    parameter is the -log10(qvalue) minimum value for calling a peak, so 2 is equivalent 
    to 0.01 or 1% FDR. The minimum peak size (length) and gap for joining nearby peaks. 
    Typically, peak size is fragment size or a little larger, and gap size is around half 
    of the fragment size.
    
        --cutoff 2 \
        --peaksize 250 \
        --peakgap 150 \
    
    Peaks are reported as [narrowPeak](http://genome.ucsc.edu/FAQ/FAQformat.html#format12)
    files, although the Macs2 `bdgpeakcall` function used in the pipeline does not 
    report pValue and qValue numbers within the file.
    
    When multiple conditions are provided, the peaks are merged into a single Bed file, 
    representing all possible peaks identified.

- Peak scoring

    The peaks are re-scored for their mean log2 fold enrichment and q-value score across
    all conditions. Additionally, normalized read counts are also collected across all 
    replicates for use in downstream differential binding analysis.
    
- Job control

	There are two options for controlling jobs and CPU usage, as some of the tools are 
	multi-threaded and many of the applications can be run concurrently. The `--job` 
	option indicates the number of simultaneous jobs that can be run concurrently. The `--cpu` 
	option indicates the number of threads available to each job. The product of the two 
	should not exceed the total number of CPUs for your machine. 

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
    
    For single-end, you can treat it as you would as ChIP-Seq.
    
- Cut site analysis

    An analysis of where the cut sites are to get a higher resolution analysis of explicitly
    open DNA (or DNase Hyper Sensitive) sites. This analysis is often chosen to 
    examine potential transcription factor binding sites that frequently occur at HS sites.
    Here, paired-end information is not required and can be treated as single-end. To 
    get a pileup of signal directly over the cut site, we provide a negative shift value 
    (to shift the alignment start in the 3' direction instead) and then extend only a 
    small distance. For example
    
        --shift=-25 \
        --extend 50 \

# Installation

HCI users running the pipeline on local servers can simply load the packages into your 
environment using a module command

    $ module load multirepchipseq
    $ multirep_macs2_pipeline.pl
    $ module unload multirepchipseq

The package only includes scripts in the `bin` directory. The scripts require additional 
library requirements of [Bio::ToolBox](https://metacpan.org/pod/Bio::ToolBox),  
[Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS), and 
[Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big), along with additional external 
applications, including [Macs2](http://github.com/taoliu/MACS/) and 
[BedTools](https://github.com/arq5x/bedtools2). Installation instructions for 
BioToolBox can be found at the [Bio::ToolBox repository](https://github.com/tjparnell/biotoolbox).
A standard [Module::Build](https://metacpan.org/pod/Module::Build) script is also 
provided for semi-automated installation of the Perl scripts.

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
