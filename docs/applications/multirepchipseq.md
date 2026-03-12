# MultiRepMacsChIPSeq - multirepchipseq

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## multirepchipseq.pl

This is a multi-threaded pipeline for processing multiple-condition,
multiple-replicate samples from ChIP-Seq, ATAC-Seq, Cut&Run, Cut&Tag, or any
other chromatin-based assays for calling peaks, particularly for differential
and/or comparative purposes. 

At least two samples (conditions) are required. Single-sample, single-replicate
analysis is not supported. If multiple replicates are available for only a
single condition, then the replicates may be treated as separate samples for
purposes of comparison.

Replicates are handled by two methods:

 1. Averaging the depth-normalized genomic coverage of the replicates into a 
    mean coverage track and calling peaks on that. This smooths the signal and
    may rescue under-performing replicates. This is always enabled.

 2. Optionally call peaks on each replicate individually, then merge the
    replicate peak calls, requiring a minimum of two overlaps to accept as a
    sample peak. This is enabled with the --indedepent option.

All sample peak calls are merged into a final call set, which are re-scored and
quantitated for comparison purposes. Multiple quality-control, comparitive, and
analytical plots are generated. An HTML report including a subset of these plots
is generated.

Multiple options are available for filtering and handling alignments, scoring,
and thresholding of peak calls. Details on these options and when to choose
them, along with numerous examples, are provided in online documentation. See
the URL below for full usage and guide.

  https://huntsmancancerinstitute.github.io/MultiRepMacsChIPSeq

Version: 21

Options:

	Input files - repeat these options for each sample (condition)
	  --chip        file1,file2...  Comma-delimited paths to experiment files
	  --name        text            Name for the sample (condition)
	  --control     file1,file2...  Optional reference control (Input), comma-delimited
	
	Output
	  --dir         directory       Directory for writing all files (./MultiRepPeakCall)
	  --out         file basename   Base filename for merged output files (all_peaks)
	
	Genome size
	  --genome      integer         Specify effective mappable genome size 
	                                (default empirically determined)
	
	Bam options
	  --pe                          Bam files are paired-end, default treat as single-end
	
	Alignment filtering options
	  --mapq        integer         Minimum mapping quality, (0)
	  --chrskip     "text"          Chromosome skip regex (chrM|MT|alt|chrUn|random|EBV|Adapter|Lambda|PhiX)
	  --blacklist   file            Coordinate (bed) file of repeats or hotspots to avoid
	                                  Default determined empirically from control samples.
	                                  Specify 'none' for no filtering.
	  --min         integer         Minimum paired-end size allowed (50 bp)
	  --max         integer         Maximum paired-end size allowed (500 bp)
	
	Duplication filtering
	  --nodedup                     Skip deduplication and use all primary, nondup alignments
	  --dupfrac     float           Target duplication rate for subsampling (0.05)
	  --maxdepth    integer         Maximum position alignment depth ()
	                                  set to 1 to remove all duplicates
	  --optdist     integer         Maximum distance for optical duplicates (0)
	                                  use 100 for HiSeq, 2500 for NovaSeq
	  --deduppair                   Run deduplication as paired-end, but coverage as single-end
	
	Fragment coverage
	  --cutsite                     Set multiple options specific for ATACSeq cutsite analysis
	  --size        integer         Predicted fragment size. REQUIRED for single-end
	  --shift       integer         Shift the fragment in special situations
	  --fraction                    Record multiple-hit alignments as fraction of hits
	  --slocal      integer         Small local lambda size (1000 bp)
	  --llocal      integer         Large local lambda size (10000 bp)
	  --cbin        integer         ChIP fragment bin size (10 bp)
	  --slbin       integer         Small local lambda bin size (50 bp)
	  --llbin       integer         Large local lambda bin size (100 bp)
	
	Chromosome-specific normalization
	  --chrnorm     float           Specific chromosome normalization factor
	  --chrapply    "text"          Apply factor to specified chromosomes via regex
	
	Peak calling
	  --independent                 Call peaks independently for each replicate and merge
	  --cutoff      number          Threshold q-value for calling peaks (2) 
	                                  Higher numbers are more significant, -1*log10(q)
	  --peaksize    integer         Minimum peak size to call (2 x fragment size)
	                                  Required for paired-end alignments.
	  --peakgap     integer         Maximum gap between peaks before merging (1 x size)
	  --broad                       Also perform broad (gapped) peak calling
	  --broadcut    number          Q-value cutoff for linking broad regions (0.5)
	  --broadgap    integer         Maximum link size between peaks in broad calls (4 x size bp)
	  --nolambda                    Skip lambda control, compare ChIP directly with control
	  --minpeakover integer         Minimum number of overlapping replicate peaks to accept
	                                  in final when merging (default n-1, minimum 2)
	  --samedepth                   Use same target depth when calculating per sample
	                                  q-value enrichment for replicate-mean peaks
	
	Peak scoring
	  --binsize     integer         Size of bins in 25 flanking peak bins for profile (100)
	  --targetdepth text            Set target depth in Millions or method for setting:
	                                  median (default), mean, min
	  --rawcounts                   Use unscaled raw counts for re-scoring peaks
	  --noplot                      Do not plot figures of results
	
	Job control
	  --cpu         integer         Number of CPUs to use per job (4)
	  --job         integer         Number of simultaneous jobs (2)
	  --dryrun                      Just print the commands without execution
	  --noorganize                  Do not organize files into subfolders when finished
	  --savebam                     Save de-duplicated bam files
	  --savebdg                     Save text bedGraph files
	
	Application  Paths
	  --bam2wig      path           (bam2wig.pl)
	  --bamdedup     path           (bam_partial_dedup.pl)
	  --bedtools     path           (bedtools)
	  --bw2bdg       path           (bigWigToBedGraph)
	  --data2wig     path           (data2wig.pl)
	  --getdata      path           (get_datasets.pl)
	  --getrel       path           (get_relative_data.pl)
	  --geteff       path           (get_chip_efficiency.pl)
	  --intersect    path           (intersect_peaks.pl)
	  --macs         path           (macs2)
	  --manwig       path           (manipulate_wig.pl)
	  --meanbdg      path           (generate_mean_bedGraph.pl)
	  --peak2bed     path           (peak2bed.pl)
	  --updatepeak   path           (update_peak_file.pl)
	  --pandoc       path           (pandoc)
	  --plotpeak     path           (plot_peak_figures.R)
	  --printchr     path           (print_chromosome_lengths.pl)
	  --reportmap    path           (report_mappable_space.pl)
	  --rscript      path           (Rscript)
	  --wig2bw       path           (wigToBigWig)



