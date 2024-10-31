# MultiRepMacsChIPSeq - recall_peaks

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

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

**NOTE**: Independently called peaks for replicates are not re-called. 
Only mean-replicate peaks are recalled. 

New peak and analysis files are placed into new subdirectories with a
numeric suffix (allowing for multiple conditions to be run
consecutively) without overwriting pre-existing files. Peaks from
different runs can subsequently be compared using the
[intersect_peaks.pl](applications/intersect_peaks.md) script, if
desired.

Version: 19

OPTIONS:

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
	 --binsize   integer           Size of bins in 25 flanking peak bins for profile (100)
	 --noplot                      Do not plot figures of results
	 
	Job control
	 --cpu       integer           Number of CPUs to use per job (4)
	 --job       integer           Number of simultaneous jobs (2)
	 --dryrun                      Just print the commands without execution
	 --noorganize                  Do not organize files into subfolders when finished
	
	Application  Paths
	 --macs        path            (macs2)
	 --manwig      path            (manipulate_wig.pl)
	 --wig2bw      path            (wigToBigWig)
	 --bw2bdg      path            (bigWigToBedGraph)
	 --printchr    path            (print_chromosome_lengths.pl)
	 --data2wig    path            (data2wig.pl)
	 --getdata     path            (get_datasets.pl)
	 --getrel      path            (get_relative_data.pl)
	 --geteff      path            (get_chip_efficiency.pl)
	 --meanbdg     path            (generate_mean_bedGraph.pl)
	 --bedtools    path            (bedtools)
	 --intersect   path            (intersect_peaks.pl)
	 --pandoc      path            (pandoc)
	 --peak2bed    path            (peak2bed.pl)
	 --updatepeak  path            (update_peak_file.pl)
	 --plotpeak    path            (plot_peak_figures.R)
	 --rscript     path            (Rscript)



