# Examples

These are example pipelines for a variety of situations. There are also 
[example alignment scripts](#example_alignment_scripts) at the end.

The associated shell scripts in this directory point to example [Bam
files](#bam_files). The Bam files are heavily sub-sampled and used here only as
working files to test that the pipeline is installed correctly and works as
expected, i.e. don't expect to interpret the results. Each script should take about
1-2 minutes to complete. Some additional parameters are included in the scripts and
excluded below for simplicity. However, they can certainly be used as templates for 
actual analysis situations.

**NOTE**: The scripts included here will not run until the project is built and/or
installed. Numerous accessory applications and R packages also need to be installed.
See the [Install](../INSTALL.md) document for details. 

To build the project and run the first example pipeline script, run the following:

    $ perl Build.PL
    $ ./Build
    $ cd examples
    $ bash ./run_se.sh


### Example ChIPSeq Pipelines


- Two sample single-end analysis

    [run_se Script](run_se.sh). Run two sample conditions, each with three replicates
    and one input. Bam files are, by default, analyzed as single-end alignments
    (irrespective of FLAG status). The mean fragment size is explicitly set to 275 bp
    (a pre-determined size). Alignments will be de-duplicated with sub-sampling to 5%
    and an optical distance of 2500 pixels. Replicates for each sample are combined
    with mean coverage and peaks called with a q-value of 0.001, minimum size
    (length) of 200 bp, and maximum gap of 100 bp. Peaks between the two samples,
    Rpd3 and Tup1, are merged into a final peak set and rescored. Comparative plots
    are generated.
    
        multirep_macs2_pipeline.pl \
        --chip Rpd3_Ch1.bam,Rpd3_Ch2.bam,Rpd3_Ch3.bam \
        --control Rpd3_Input.bam \
        --name Rpd3 \
        --chip Tup1_Ch1.bam,Tup1_Ch2.bam,Tup1_Ch3.bam \
        --control Tup1_Input.bam \
        --name Tup1 \
        --size 275 \
        --out se_all \
        --dupfrac 0.05 \
        --optdist 2500 \
        --cutoff 3 \
        --peaksize 200 \
        --peakgap 100 


- Two sample paired-end analysis

    [run_pe Script](run_pe.sh). Same as above, but indicate that the Bam files are
    paired-end alignments by adding the `--pe` option. Only properly-paired
    alignments are analyzed; singletons and non-paired alignments are silently
    skipped. The fragment size is not necessary (unless the `--peaksize` option was
    not defined).
    
        multirep_macs2_pipeline.pl \
        --chip Rpd3_Ch1.bam,Rpd3_Ch2.bam,Rpd3_Ch3.bam \
        --control Rpd3_Input.bam \
        --name Rpd3 \
        --chip Tup1_Ch1.bam,Tup1_Ch2.bam,Tup1_Ch3.bam \
        --control Tup1_Input.bam \
        --name Tup1 \
        --pe \
        --out pe_all \
        --dupfrac 0.05 \
        --optdist 2500 \
        --cutoff 3 \
        --peaksize 200 \
        --peakgap 100 

- One sample paired-end with individual replicates 

    [run_individual_replicates Script](run_individual_replicates.sh). This has only
    one ChIP sample with three replicates. The replicates are indicated separately,
    which will allow peaks to be called individually, merged, and comparative plots
    generated.

        multirep_macs2_pipeline.pl \
        --chip Tup1_Ch1.bam \
        --name Tup1_Ch1 \
        --chip Tup1_Ch2.bam \
        --name Tup1_Ch2 \
        --chip Tup1_Ch3.bam \
        --name Tup1_Ch3 \
        --control Tup1_Input.bam \
        --pe \
        --out tup1_all \
        --cutoff 3 \
        --peaksize 200 \
        --peakgap 100 


- Two Sample paired-end with independent replicates

    [run_independent Script](run_independent_pe.sh). This example is with two
    paired-end samples but with unequal numbers of replicates. The `--independent`
    flag is used to call peaks independently for each replicate, then merged into a
    sample-level peak file. Note that, by default, a minimum of two overlapping
    replicate peaks will be required for Rpd3 since three are provided. In addition
    to replicate-merged peaks, mean-replicate peaks will also be called by default.
    Plots generated for each analysis are placed into separate folders. 

        multirep_macs2_pipeline.pl \
        --chip Rpd3_Ch1.bam,Rpd3_Ch2.bam,Rpd3_Ch3.bam \
        --control Rpd3_Input.bam \
        --name Rpd3 \
        --chip Tup1_Ch1.bam \
        --control Tup1_Input.bam \
        --name Tup1 \
        --pe \
        --out pe_all \
        --independent \
        --cutoff 3 \
        --peaksize 200 \
        --peakgap 100 \


- Two Sample paired-end with independent replicates and no deduplication

    [run_independent_nodup Script](run_independent_nodup.sh). This example is with two
    ChIP samples, each with three replicates, and with de-duplication turned off.
    De-duplication should be explicitly turned off when the bam files have already 
    been de-duplicated, particularly when using 
    [UMI-based deduplication](https://github.com/HuntsmanCancerInstitute/UMIScripts)
    which leaves behind duplicates of biological origin. Even without de-duplication,
    Bam files are automatically filtered prior to independent peak calling. 

        multirep_macs2_pipeline.pl \
        --chip Rpd3_Ch1.bam,Rpd3_Ch2.bam,Rpd3_Ch3.bam \
        --control Rpd3_Input.bam \
        --name Rpd3 \
        --chip Tup1_Ch1.bam,Tup1_Ch2.bam,Tup1_Ch3.bam \
        --control Tup1_Input.bam \
        --name Tup1 \
        --pe \
        --out all_nodup \
        --nodedup \
        --independent \
        --cutoff 3 \
        --peaksize 200 \
        --peakgap 100 

- Two Sample Paired-end without controls

    [run_pe_nocontrol Script](run_pe_nocontrol.sh). This example is with two
    paired-end samples, each with three replicates, but without reference controls. A
    global mean coverage is generated as the _lambda_ control from the ChIP samples.

        multirep_macs2_pipeline.pl \
        --chip Rpd3_Ch1.bam,Rpd3_Ch2.bam,Rpd3_Ch3.bam \
        --name Rpd3 \
        --chip Tup1_Ch1.bam,Tup1_Ch2.bam,Tup1_Ch3.bam \
        --name Tup1 \
        --pe \
        --out all_nocontrol \
        --cutoff 3 \
        --peaksize 200 \
        --peakgap 100 


- Two Sample Paired-end broad call

    [run_pe_gapped Script](run_pe_gapped.sh). This example is with two paired-end
    samples with broad (or gapped) peak calling. The tolerated gap is set to 1000 bp
    and the minimum q-value tolerated is 1, or 0.1. Narrow peak calling is always
    performed. 

        multirep_macs2_pipeline.pl \
        --chip Rpd3_Ch1.bam,Rpd3_Ch2.bam,Rpd3_Ch3.bam \
        --control Rpd3_Input.bam \
        --name Rpd3 \
        --chip Tup1_Ch1.bam,Tup1_Ch2.bam,Tup1_Ch3.bam \
        --control Tup1_Input.bam \
        --name Tup1 \
        --pe \
        --out pe_all \
        --cutoff 3 \
        --peaksize 200 \
        --peakgap 100 \
        --broad \
        --broadgap 1000 \
        --broadcut 1 \


- ATAC-Seq cut site analysis

    [run_cutsite Script](run_cutsite.sh). This example runs a cut site analysis
    typical of ATAC-Seq paired-end samples. Alignments are de-duplicated as
    paired-end (indicated by the `--deduppair` option), but analysis is performed as
    single-end, with 100 bp coverage centered over the 5' end by shifting the
    coverage -50 bp. Additionally, peak sizes are adjusted to account for narrower
    peaks than ChIP.

        multirep_macs2_pipeline.pl \
        --chip Rpd3_Ch1.bam,Rpd3_Ch2.bam,Rpd3_Ch3.bam \
        --name Rpd3 \
        --chip Tup1_Ch1.bam,Tup1_Ch2.bam,Tup1_Ch3.bam \
        --name Tup1 \
        --out all_cut \
        --deduppair \
        --dupfrac 0.1 \
        --shift -50 \
        --size 100 \
        --cutoff 2 \
        --peaksize 150 \
        --peakgap 50 


- Single-end and Paired-end combination analysis

    [run_se_pe_merge_recall Script](run_se_pe_merge_recall.sh). This script example
    runs two pipelines sequentially, one sample as paired-end and the other sample as
    single-end. The resulting files are then merged into a third folder to mimic the
    samples having been run in a single run. Finally, peaks are recalled using the
    [recall_peaks](../Applications.md#recall_peakspl) script to generate a merged
    peak comparison analysis. While this may not be an ideal situation, it is an 
    illustration of what could be achieved in weird situations.




### Bam files

The example Bam files are derived from the publication 
[Parnell et al, 2020, PLoS Genetics](https://pubmed.ncbi.nlm.nih.gov/33382702/).
The raw Fastq data is archived in GEO Accession 
[GSE158180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158180).

These are two ChIPSeq experiments from _Saccharomyces cerevisiae_ of two different
chromatin factors, a corepressor (Tup1) and a histone deacetylase (Rpd3). There are
three biological replicates and one Input reference per experiment.

|FileName       | SampleName    | GEO ID    |
|---------------|---------------|-----------|
|Tup1_Ch1.bam   | Tup1-V5_rep1  | GSM4794781|
|Tup1_Ch2.bam   | Tup1-V5_rep2  | GSM4794782|
|Tup1_Ch3.bam   | Tup1-V5_rep3  | GSM4794783|
|Rpd3_Ch1.bam   | Rpd3-V5_rep1  | GSM4794784|
|Rpd3_Ch2.bam   | Rpd3-V5_rep2  | GSM4794785|
|Rpd3_Ch3.bam   | Rpd3-V5_rep3  | GSM4794786|
|Tup1_Input.bam | Tup1-V5_Input | GSM4794788|
|Rpd3_Input.bam | Rpd3-V5_Input | GSM4794789|

Example Bam files were aligned to SacCer3. Properly paired alignments on chrI were 
extracted and subsampled to < 5%. Read sequence, qualities, and CIGAR string were 
condensed to the first aligned base. Read names were condensed to tile coordinates 
only. Extraneous SAM attributes were excluded.



### Example Alignment Scripts

These are example alignment script templates for aligning raw Fastq files to generate 
alignment Bam files suitable for use in the pipeline. They utilize 
[Novocraft novoAlign](https://www.novocraft.com/products/novoalign/) for alignment, 
and run a few metrics applications to determine duplication rate and fragment size. 

For users of the HCI 
[pysano](https://uofuhealth.utah.edu/huntsman/shared-resources/gba/bioinformatics/infrastructure/pysano.php) 
system for executing jobs at [CHPC](https://www.chpc.utah.edu), you can turn these into 
a `cmd.txt` file by simply removing the first three lines. Applications from packages 
installed under a modules environment are called using a simple wrapper script in the 
HCI production environment, e.g. `$APP/modulesoftware/bam_partial_dedup`. 

- Single-end alignment

    [align_singleEnd_example Script](align_singleEnd_example.sh). An example script 
    to align single-end Fastq files and determine fragment size using both
    Macs2 `predictd` function and BioToolBox 
    [bam2wig](https://metacpan.org/pod/bam2wig.pl) as described. Duplication rate 
    is determined using [bam_partial_dedup.pl](../applications.md#bam_partial_deduppl).

- Paired-end alignment

    [align_pairedEnd_example Script](align_pairedEnd_example.sh). An example 
    script to align paired-end Fastq files, either ChIPSeq or ATAC (Tn5 digested) 
    samples. Duplication rate is determined using 
    [bam_partial_dedup.pl](../applications.md#bam_partial_deduppl). Insertion size 
    is checked with 
    [Picard CollectInsertSizeMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics).
    
- Paired-end UMI alignment

    [align_pairedEnd_UMI_example Script](align_pairedEnd_UMI_example.sh). This 
    example script is for samples prepared with a library with Unique Molecular 
    Indexes (UMIs) to uniquely identify PCR duplicates. The UMI codes are provided 
    in a third Fastq and are merged as SAM comments into new Fastq files using the 
    `merge_umi_fastq.pl` script from the 
    [UMIScripts package](https://github.com/HuntsmanCancerInstitute/UMIScripts).
    Duplicate alignments are then removed using the `bam_umi_dedup.pl` script 
    from the same package (alignments may also be simply marked without removal 
    using the `--mark` option). Insertion size is checked with 
    [Picard CollectInsertSizeMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics).
    

You can combine the metrics with the script 
[combine_replicate_data.pl](../applications.md#combine_replicate_datapl). It will 
look for and read `stdout.txt` and `stderr.txt` files.


