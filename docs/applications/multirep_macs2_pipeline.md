# MultiRepMacsChIPSeq - multirep\_macs2\_pipeline

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## multirep\_macs2\_pipeline.pl

This is a wrapper for calling and/or comparing peaks in ChIPSeq or ATACSeq with
single or multiple replicas using the [Macs2](https://pypi.org/project/MACS2/)
ChIPSeq caller. It uses [BioToolBox](http://tjparnell.github.io/biotoolbox)
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

Version: 20

Options:

	Input files
	  --chip        file1,file2...  Repeat for each sample set
	  --name        text            Repeat for each sample
	  --control     file1,file2...  Repeat if matching multiple samples
	
	Output
	  --dir         directory       Directory for writing all files (./MultiRepPeakCall)
	  --out         file basename   Base filename for merged output files (merged)
	
	Genome size
	  --genome      integer         Specify effective mappable genome size 
	                                (default empirically determined)
	
	Bam options
	  --pe                          Bam files are paired-end, default treat as single-end
	
	Alignment filtering options
	  --mapq        integer         Minimum mapping quality, (0)
	  --chrskip     "text"          Chromosome skip regex (chrM|MT|alt|Adapter|Lambda|PhiX)
	  --blacklist   file            Bed file of repeats or hotspots to avoid
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
	  --targetdepth text            Set method for sequence depth scaling for all count data:
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



