#!/usr/bin/perl

# Timothy J. Parnell, PhD
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112
#  
# This package is free software; you can redistribute it and/or modify
# it under the terms of the Artistic License 2.0.  
# 
# Updated versions of this file may be found in the repository
# https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq


use strict;
use IO::File;
use File::Copy;
use File::Spec;
use File::Which;
use File::Path qw(make_path);
use List::Util qw(min uniqstr);
use Getopt::Long;
use Parallel::ForkManager;
use Bio::ToolBox::utility qw(simplify_dataset_name);

my $VERSION = 15.2;

# parameters
my %opts = (
	dir         => './',
	out         => 'merged',
	genome      => 0,
	species     => '',
	paired      => 0,
	mapq        => 0,
	fraction    => 0,
	minsize     => 50,
	maxsize     => 500,
	dedup       => 1,
	maxdup      => undef, # old option
	maxdepth    => undef,
	dupfrac     => 0.05,
	optdist     => 0,
	deduppair   => 0,
	savebam     => 0,
	fragsize    => 250,
	shiftsize   => 0,
	slocal      => 1000,
	llocal      => 10000,
	cutoff      => 2,
	targetdep   => undef,
	peaksize    => undef,
	peakgap     => undef,
	broad       => 0,
	broadcut    => 0.5,
	broadgap    => undef,
	linkqv      => 1,
	gaplink     => undef,
	lambda      => 1,
	chrskip     => "chrM|MT|lambda|Adapter|PhiX",
	blacklist   => undef,
	cpu         => 4,
	job         => 2,
	chipbin     => 10,
	slocalbin   => 50,
	llocalbin   => 100,
	chrapply    => undef,
	rawcounts   => 0,
	savebdg     => 0,
	binsize     => 40,
	genomewin   => 0,
	discard     => 10,
	repmean     => 0,
	plot        => 1,
	dryrun      => 0,
	organize    => 1,
	bam2wig     => sprintf("%s", which 'bam2wig.pl'),
	bamdedup    => sprintf("%s", which 'bam_partial_dedup.pl'),
	macs        => sprintf("%s", which 'macs2'),
	manwig      => sprintf("%s", which 'manipulate_wig.pl'),
	wig2bw      => sprintf("%s", which 'wigToBigWig'),
	bw2bdg      => sprintf("%s", which 'bigWigToBedGraph'),
	bedtools    => sprintf("%s", which 'bedtools'),
	getdata     => sprintf("%s", which 'get_datasets.pl'),
	getrel      => sprintf("%s", which 'get_relative_data.pl'),
	geteff      => sprintf("%s", which 'get_chip_efficiency.pl'),
	printchr    => sprintf("%s", which 'print_chromosome_lengths.pl'),
	data2wig    => sprintf("%s", which 'data2wig.pl'),
	meanbdg     => sprintf("%s", which 'generate_mean_bedGraph.pl'),
	intersect   => sprintf("%s", which 'intersect_peaks.pl'),
	peak2bed    => sprintf("%s", which 'peak2bed.pl'),
	combrep     => sprintf("%s", which 'combine_replicate_data.pl'),
	plotpeak    => sprintf("%s", which 'plot_peak_figures.R'),
	rscript     => sprintf("%s", which 'Rscript'),
	reportmap   => sprintf("%s", which 'report_mappable_space.pl'),
);
my @names;
my @chips;
my @controls;
my @chip_scales;
my @control_scales;
my @chrnorms;

my $help = 0;
my $documentation = <<DOC;

======= ChIP Wrapper ==========

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
fragment (or d in Macs2 parlance), small lambda (default $opts{slocal} bp), and 
large lambda (default $opts{llocal} bp) fragment coverage. If desired, either small or 
local lambda may be turned off by setting to 0. To completely turn off lambda, set the 
nolambda option, whereupon only the control fragment is directly used as reference. 
If no control file is provided, then the global mean from the ChIP file is used 
as a (poor) substitute. 

Advanced users may provide one processed bigWig file per ChIP or control sample. 

See https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq for full 
usage and guide.

Version: $VERSION

Options:
 Input files
  --chip      file1,file2...    Repeat for each sample set
  --name      text              Repeat for each sample
  --control   file1,file2...    Repeat if matching multiple samples
 
 Output
  --dir       directory         Directory for writing all files ($opts{dir})
  --out       file basename     Base filename for merged output files ($opts{out})
  
 Genome size
  --genome    integer           Specify effective mappable genome size 
                                  (default empirically determined)
  
 Bam options
  --mapq      integer           Minimum mapping quality, ($opts{mapq})
  --pe                          Bam files are paired-end, default treat as single-end
  --min       integer           Minimum paired-end size allowed ($opts{minsize} bp)
  --max       integer           Maximum paired-end size allowed ($opts{maxsize} bp)
  --fraction                    Record multiple-hit alignments as fraction of hits
 
 Bam filtering options
  --chrskip   "text"            Chromosome skip regex ($opts{chrskip})
  --blacklist file              Bed file of repeats or hotspots to avoid
                                  Determined empirically from control (Input) samples
  
 Duplication filtering
  --nodedup                     Skip deduplication and take everything as is
  --dupfrac   float             Target duplication rate for subsampling ($opts{dupfrac})
  --maxdepth  integer           Maximum position alignment depth ($opts{maxdepth})
                                  set to 1 to remove all duplicates
  --optdist   integer           Maximum distance for optical duplicates ($opts{optdist})
                                  use 100 for HiSeq, 2500 for NovaSeq
  --deduppair                   Run deduplication as paired-end, but coverage as single-end
                                  e.g. for ATAC-Seq cut site analysis
  --savebam                     Save de-duplicated bam files

 Fragment coverage
  --size      integer           Predicted fragment size. REQUIRED for single-end
  --shift     integer           Shift the fragment, e.g. ATACSeq ($opts{shiftsize} bp)
  --slocal    integer           Small local lambda size ($opts{slocal} bp)
  --llocal    integer           Large local lambda size ($opts{llocal} bp)
  --cbin      integer           ChIP fragment bin size ($opts{chipbin} bp)
  --slbin     integer           Small local lambda bin size ($opts{slocalbin} bp)
  --llbin     integer           Large local lambda bin size ($opts{llocalbin} bp)

 Chromosome-specific normalization
  --chrnorm   float             Specific chromosome normalization factor
  --chrapply  "text"            Apply factor to specified chromosomes via regex
 
 Peak calling
  --cutoff    number            Threshold q-value for calling peaks ($opts{cutoff}) 
                                  Higher numbers are more significant, -1*log10(q)
  --peaksize  integer           Minimum peak size to call (2 x fragment size)
                                  Required for paired-end alignments.
  --peakgap   integer           Maximum gap between peaks before merging (1 x size)
  --broad                       Also perform broad (gapped) peak calling
  --broadcut  number            Q-value cutoff for linking broad regions ($opts{broadcut})
  --broadgap  integer           Maximum link size between peaks in broad calls (4 x size bp)
  --nolambda                    Skip lambda control, compare ChIP directly with control
  --savebdg                     Save q-value bdg files for further custom calling
  
 Peak scoring
  --binsize   integer           Size of bins in 25 flanking peak bins for profile ($opts{binsize})
  --window    integer           Collect counts across genome in given window size
  --discard   number            Discard genome windows with replicate sum below number ($opts{discard})
  --rawcounts                   Use unscaled raw counts for re-scoring peaks
  --repmean                     Combine replicate counts as mean for each sample set
  --noplot                      Do not plot figures of results
  
 Job control
  --cpu       integer           Number of CPUs to use per job ($opts{cpu})
  --job       integer           Number of simultaneous jobs ($opts{job})
  --dryrun                      Just print the commands without execution
  --noorganize                  Do not organize files into subfolders when finished

 Application  Paths
  --bam2wig   path             ($opts{bam2wig})
  --bamdedup  path             ($opts{bamdedup})
  --macs      path             ($opts{macs})
  --manwig    path             ($opts{manwig})
  --wig2bw    path             ($opts{wig2bw})
  --bw2bdg    path             ($opts{bw2bdg})
  --printchr  path             ($opts{printchr})
  --data2wig  path             ($opts{data2wig})
  --getdata   path             ($opts{getdata})
  --getrel    path             ($opts{getrel})
  --geteff    path             ($opts{geteff})
  --meanbdg   path             ($opts{meanbdg})
  --bedtools  path             ($opts{bedtools})
  --intersect path             ($opts{intersect})
  --peak2bed  path             ($opts{peak2bed})
  --combrep   path             ($opts{combrep})
  --plotpeak  path             ($opts{plotpeak})
  --rscript   path             ($opts{rscript})
  --reportmap path             ($opts{reportmap})
DOC


### Inputs
unless (@ARGV) {
	print <<HELP;

======= ChIP Wrapper ==========

This is a wrapper for calling and/or comparing peaks in ChIPSeq or ATACSeq with single 
or multiple replicas using the Macs2 ChIPSeq caller. It uses BioToolBox applications to 
normalize duplicate levels and read depths between samples and replicates.

See https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq for full 
usage and guide.

Use --help to display options.

Version: $VERSION

HELP
	exit 1;
}

GetOptions(
	'dir=s'                 => \$opts{dir},
	'out=s'                 => \$opts{out},
	'chip=s'				=> \@chips,
	'control=s'             => \@controls,
	'name=s'                => \@names,
	'chscale=s'             => \@chip_scales,
	'coscale=s'             => \@control_scales,
	'species=s'             => \$opts{species},
	'genome=i'              => \$opts{genome},
	'mapq=i'                => \$opts{mapq},
	'paired|pe!'            => \$opts{paired},
	'min=i'                 => \$opts{minsize},
	'max=i'                 => \$opts{maxsize},
	'fraction!'             => \$opts{fraction},
	'chrskip=s'             => \$opts{chrskip},
	'blacklist=s'           => \$opts{blacklist},
	'dedup!'                => \$opts{dedup},
	'dupfrac=f'             => \$opts{dupfrac},
	'maxdup=i'              => \$opts{maxdup}, # old option
	'maxdepth=i'            => \$opts{maxdepth},
	'optdist=i'             => \$opts{optdist},
	'deduppair!'            => \$opts{deduppair},
	'savebam!'              => \$opts{savebam},
	'size=i'                => \$opts{fragsize},
	'shift=i'               => \$opts{shiftsize},
	'slocal=i'              => \$opts{slocal},
	'llocal=i'              => \$opts{llocal},
	'cbin=i'                => \$opts{chipbin},
	'slbin=i'               => \$opts{slocalbin},
	'llbin=i'               => \$opts{llocalbin},
	'chrnorm=f'             => \@chrnorms,
	'chrapply=s'            => \$opts{chrapply},
	'cutoff=f'              => \$opts{cutoff},
	'tdep=f'                => \$opts{targetdep},
	'peaksize=i'            => \$opts{peaksize},
	'peakgap=i'             => \$opts{peakgap},
	'broad!'                => \$opts{broad},
	'broadcut=f'            => \$opts{broadcut},
	'broadgap=i'            => \$opts{broadgap},
	'lambda!'               => \$opts{lambda},
	'savebdg!'              => \$opts{savebdg},
	'window=i'              => \$opts{window},
	'discard=f'             => \$opts{discard},
	'repmean!'              => \$opts{repmean},
	'plot!'                 => \$opts{plot},
	'cpu=i'                 => \$opts{cpu},
	'job=i'                 => \$opts{job},
	'dryrun!'               => \$opts{dryrun},
	'organize!'             => \$opts{organize},
	'help!'                 => \$help,
	'bam2wig=s'             => \$opts{bam2wig},
	'bamdedup=s'            => \$opts{bamdedup},
	'macs=s'                => \$opts{macs},
	'manwig=s'              => \$opts{manwig},
	'wig2bw=s'              => \$opts{wig2bw},
	'bw2bdg=s'              => \$opts{bw2bdg},
	'printchr=s'            => \$opts{printchr},
	'data2wig=s'            => \$opts{data2wig},
	'getdata=s'             => \$opts{getdata},
	'getrel=s'              => \$opts{getrel},
	'geteff=s'              => \$opts{geteff},
	'meanbdg=s'             => \$opts{meanbdg},
	'bedtools=s'            => \$opts{bedtools},
	'intersect=s'           => \$opts{intersect},
	'peak2bed=s'            => \$opts{peak2bed},
	'combrep=s'             => \$opts{combrep},
	'plotpeak=s'            => \$opts{plotpeak},
	'rscript=s'             => \$opts{rscript},
	'reportmap=s'           => \$opts{reportmap},
) or die "unrecognized option(s)!\n";

if ($help) {
	print $documentation;
	exit 1;
}



### Begin main pipeline
print "======== ChIPSeq multi-replicate pipeline ==========\n";
print "\nversion $VERSION\n";
my $start = time;



# steps
my @finished_commands;
check_inputs();
print_start();
check_input_files();
my @Jobs = generate_job_file_structure();
my $progress_file = File::Spec->catfile($opts{dir}, $opts{out} . '.progress.txt');
my %progress = check_progress_file();
my $chromofile = generate_chr_file();
run_input_peak_detection();
run_dedup();
run_bam_check();
run_mappable_space_report();
run_bam_fragment_conversion();
run_bam_count_conversion();
run_lambda_control();
run_bw_conversion();
run_bdgcmp();
run_call_peaks();
run_clean_peaks();
run_peak_merge();
run_bdg_conversion();
run_rescore();
run_efficiency();
run_plot_peaks();
run_cleanup();
run_organize();

# final statement
printf "\n\nFinished in %.1f minutes\n", (time -$start) / 60;
print "======== Finished ChIPSeq multi-replicate pipeline ==========\n";




############### Subroutines ########################################################

sub check_inputs {
	if (@ARGV) {
		die sprintf("There are unrecognized leftover items on the command line!\n Did you leave spaces in your --chip or --control file lists?\nItems:\n %s\n", join("\n ", @ARGV));
	}
	unless (@chips) {
		die "No ChIP file(s) defined!\n";
	}
	unless (@names) {
		die "No name(s) defined!\n";
	}
	unless (scalar(@chips) == scalar(@names)) {
		die "Unequal number of ChIP samples and names!\n";
	}
	if (scalar(@controls) > 1 and scalar(@controls) != scalar(@chips)) {
		die "Unequal number of control and ChIP samples!\n";
	}
	elsif (scalar(@controls) == 0) {
		# no controls, turn off lambda
		$opts{slocal} = 0;
		$opts{llocal} = 0;
		$opts{lambda} = 0;
	}
	if (scalar(@chip_scales) or scalar(@control_scales)) {
		# no longer recommended
			print <<MESSAGE;

WARNING: Manually setting ChIP and/or Control scaling factors!!!!
Please be aware that manually scaling coverage depth is an advanced option
and should only be done when you're aware of the ramifications. Inappropriate
scaling may artificially alter the expected statistics required for accurate
peak calling and may reduce confidence of identified peaks.

MESSAGE
	}
	if (scalar(@chip_scales) and scalar(@chip_scales) != scalar(@chips)) {
		die "Unequal number of ChIP samples and ChIP scale factors!\n";
	}
	if (scalar(@control_scales) and scalar(@control_scales) != scalar(@controls)) {
		die "Unequal number of control samples and control scale factors!\n";
	}
	if (scalar(@chrnorms) and not $opts{chrapply}) {
		die "Chromosome normalization factors given but no chromosome specified!\n";
	}
	if (not scalar(@chrnorms) and $opts{chrapply}) {
		die "Chromosome name for normalization specified but no factors given!\n";
	}
	if (scalar(@chrnorms) and scalar(@chrnorms) != scalar(@chips)) {
		# apply to all the ChIPs
		if (scalar @chrnorms == 1) {
			print "Using the same chromosome normalization factor for each ChIP sample\n";
			my $n = shift @chrnorms;
			foreach (@names) {
				push @chrnorms, $n;
			}
		}
		else {
			die "Unequal number of chromosome normalization factors and ChIP samples!\n";
		}
	}
	if (scalar(@chrnorms) > 1 and scalar(@controls) == 1) {
		print "Using first chromosome normalization factor for universal control!\n";
	}
	if ($opts{species}) {
		print <<MESSAGE;

WARNING: Specifiying species is now deprecated. The genome mappable size is 
now determined empirically from all provided Bam files using the script 
report_mappable_space.pl, or an explicit genome mappable size may be provided 
with the --genome option. 

MESSAGE
	}
	
	# directory
	unless ($opts{dryrun}) {
		unless (-e $opts{dir} and -d _) {
			make_path($opts{dir}) or die sprintf("unable to make directory %s! $!\n", $opts{dir});
		}
		unless (-w $opts{dir}) {
			die sprintf("target directory %s is not writable!\n", $opts{dir});
		}
	}
	# target depth
	if (defined $opts{targetdep}) {
		my $files = join(",", @chips, @controls);
		if ($files =~ /\.bam/i) {
			# no longer recommend manually setting with bam files
			# bigWig files would be acceptable
			print <<MESSAGE;

WARNING: Manually setting the target sequence depth is not recommended!!!!
Target depth is now automatically calculated empirically as the minimum
number of accepted fragments from all provided bam files. Reconsider using
this option unless you understand the ramifications.

MESSAGE
		}
	}
	
	# sizes
	if ($opts{fragsize}) {
		$opts{peaksize} ||= 2 * $opts{fragsize};
	}
	else {
		# no fragment size defined? might be ok
		if ($opts{paired}) {
			# paired fragments
			if (not $opts{peaksize}) {
				print "\n WARNING! Setting minimum peak size to 500 bp, but this should be manually\n set based on mean alignment insert size and nature of experiment.\n\n";
				$opts{peaksize} = 500;
				$opts{fragsize} = 250; # for expected background normalization
			}
		}
		else {
			die " Must set an estimated mean fragment size for single-end alignments!\n  Run 'macs2 predictd' or 'bam2wig.pl --shift --model'\n";
		}
	}
	unless (defined $opts{peakgap}) {
		$opts{peakgap} = $opts{fragsize} ? $opts{fragsize} : int($opts{peaksize}/2);
	}
	unless (defined $opts{broadgap}) {
		$opts{broadgap} = $opts{fragsize} ? 4 * $opts{fragsize} : 2 * $opts{peaksize};
	}
	
	# add parameters to option hash for printing configuration
	$opts{chipscale} = join(", ", @chip_scales);
	$opts{controlscale} = join(", ", @control_scales);
	$opts{chrnorms} = join(", ", @chrnorms);
	
	# exclusion list
	if ($opts{blacklist} and not -e $opts{blacklist}) {
		printf("\n WARNING! Unable to find specified black list file '%s'!\n", $opts{blacklist});
		if (scalar(@controls)) {
			print " Defaulting to using input-derived exclusion list\n";
			$opts{blacklist} = 'input';
		}
		else {
			undef $opts{blacklist};
		}
	}
	if (not defined $opts{blacklist} and scalar(@controls)) {
		$opts{blacklist} = 'input';
	}
	
	# max depth-duplication confusion
	if (defined $opts{maxdup}) {
		# because this was inappropriately named before
		print " \n WARNING: The --maxdup option is now --maxdepth\n";
		$opts{maxdepth} = $opts{maxdup};
	}
	delete $opts{maxdup};
}

sub print_start {
	if ($opts{dryrun}) {
		print <<DRYRUN;

======= Dry Run =======
Some values are empirically determined during execution and are made up here.
Some files may not be generated during actual execution, but commands should
mostly be complete. No files will be generated.
=======================

DRYRUN
	}
	print "\n\n======= Samples\n";
	for my $i (0 .. $#names) {
		printf " %s: %s\n", $names[$i], $chips[$i]; 
		printf " Control: %s\n", $controls[$i] if $i < scalar(@controls);
	}
	print "\n\n======= Configuration\n";
	foreach my $k (sort {$a cmp $b} keys %opts) {
		printf "%12s  %s\n", $k, $opts{$k};
	}
	print "\n\n";
}

sub check_input_files {
	print "\n\n======= Checking input files\n";
	my $error = 0;
	foreach my $i (0 .. $#names) {
		my @list = split ',', $chips[$i];
		if (scalar(@list) != uniqstr(@list)) {
			printf " Duplicate entries found in %s ChIP entries!\n", $names[$i];
			# fix it
			@list = uniqstr(@list);
			$chips[$i] = join(',', @list);
		}
		foreach my $f (@list) {
			unless (-e $f) {
				printf(" can't find %s ChIP file %s!\n", $names[$i], $f);
				$error++;
			}
		}
		if ($controls[$i]) {
			my @list2 = split ',', $controls[$i];
			if (scalar(@list2) != uniqstr(@list2)) {
				printf " Duplicate entries found in %s control entries!\n", $names[$i];
				# fix it
				@list2 = uniqstr(@list2);
				$chips[$i] = join(',', @list2);
			}
			foreach my $f (@list2) {
				unless (-e $f) {
					printf(" can't find %s control file %s!\n", $names[$i], $f);
					$error++;
				}
			}
		}
	}
	if ($error) {
		die "ERROR! Missing $error input files!\n";
	}
	else {
		print " All input files found\n";
	}
}

sub generate_job_file_structure {
	# we pass this off to ChIPjob package with 6 options: 
	# name, chip files, control files, chip scale factor, control scale factor, 
	# chromosome normalization factor
	my @jobs;
	
	# first check for a unversal control
	if (scalar(@controls) == 1 and scalar(@chips) > 1) {
		print "Using only one control for multiple ChIP experiments\n";
		my $universal_control = shift @controls;
		
		# generate universal name
		my $universal_name;
		if ($universal_control =~ /\.bam$/i) {
			# one bam file
			$universal_name = $opts{out} . '_control';
			if ($opts{lambda}) {
				# add lambda control if we're using that, otherwise leave it be
				$universal_name .= '.lambda_control';
			}
			else {
				$universal_name .= '.control_fragment';
			}
		}
		elsif ($universal_control =~ /\.(?:bw|bigwig)$/i) {
			# a pre-processed bigWig file
			(undef, undef, $universal_name) = File::Spec->splitpath($universal_control);
			$universal_name =~ s/\.(?:bw|bigwig)$//;
		}
		else {
			die "unrecognized control file '$universal_control'!\n";
		}
		my $universal_scale = shift @control_scales || undef;
		
		# add back universal control with special prefix for each ChIP job
		foreach (@names) {
			push @controls, 'Custom-Universal-' . $universal_name;
		}
		
		# generate job
		push @jobs, ChIPjob->new($universal_name, '', $universal_control, undef, 
			$universal_scale, $chrnorms[0] || undef);
	}
	
	# walk through each given job
	for my $i (0 .. $#names) {
		push @jobs, ChIPjob->new($names[$i], $chips[$i], $controls[$i] || undef, 
			$chip_scales[$i] || undef, $control_scales[$i] || undef, 
			$chrnorms[$i] || undef );
	}
	
	return @jobs;
}

sub check_progress_file {
	my %p = (
		control_peak  => 0,
		deduplication => 0,
		mappable_size => 0,
		fragment      => 0,
		count         => 0,
		lambda        => 0,
		bw2bdg        => 0,
		bdgcmp        => 0,
		callpeak      => 0,
		cleanpeak     => 0,
		bdg2bw        => 0,
		peakmerge     => 0,
		rescore       => 0,
		efficiency    => 0
	);
	if (-e $progress_file) {
		my $fh = IO::File->new($progress_file, '<');
		while (my $line = $fh->getline) {
			chomp $line;
			$p{$line} = 1 if exists $p{$line};
			if ($line eq 'bam2wig') {
				# old value!? shouldn't happen but just in case
				$p{fragment} = 1;
				$p{count} = 1;
			}
		}
		$fh->close;
	}
	return %p;
}

sub generate_chr_file {
	my @bams = split(',', $chips[0]);
	my $example = shift @bams;
		# this will work regardless if example is bam or bigWig
	my $chromofile = File::Spec->catfile($opts{dir},"chrom_sizes.temp.txt");
	return $chromofile if $opts{dryrun}; # just pretend
	if (-e $chromofile) {
		return $chromofile;
	}
	unless ($opts{printchr}) {
		die "no print_chromosome_lengths.pl script in path!\n";
	}
	system(sprintf("%s --db %s --chrskip '%s' --out %s", 
		$opts{printchr}, $example, $opts{chrskip}, $chromofile));
	die "chromosome file $chromofile could not be generated!\n" unless -e $chromofile;
	return $chromofile;
}

sub execute_commands {
	my $commands = shift;
	printf "Excecuting %d commands\n", scalar @$commands;
	
	# dry run
	if ($opts{dryrun}) {
		# we just go through the motions here
		foreach my $command (@$commands) {
			printf "=== Job: %s\n", $command->[0];
		}
		push @finished_commands, @$commands;
		return;
	}
	
	# execute jobs
	if ($opts{job} > 1) {
		my $pm = Parallel::ForkManager->new($opts{job});
		foreach my $command (@$commands) {
			next if check_command_finished($command, 1);
			printf "=== Job: %s\n", $command->[0];
			
			# check for simple rm commands
			if ($command->[0] =~ /^rm (.+)$/) {
				# we don't need to fork a new process just to execute a rm command
				unlink($1);
				next;
			}
			
			# fork to execute
			$pm->start and next;
			# in child
			system($command->[0]);
			$pm->finish;
		}
		$pm->wait_all_children;
	}
	else {
		foreach my $command (@$commands) {
			next if check_command_finished($command, 1);
			printf "=== Job: %s\n", $command->[0];
			system($command->[0]);
		}
	}
	
	# check that commands actually produced something
	my @errors;
	foreach my $command (@$commands) {
		# returns a non-zero function if command didn't run
		push @errors, $command unless check_command_finished($command, 0);
	}
	if (@errors) {
		print "\n\n ======= Errors ======\n";
		print " The following jobs did not generate expected output\n";
		foreach (@errors) {
			printf "=== ERROR: %s\n", $_->[0];
		}
		die "\nCheck log files for errors\n";
	}
	
	push @finished_commands, @$commands;
}

sub check_command_finished {
	my ($command, $talk) = @_;
	# returns true if command appears finished
	
	# command bits
	my ($command_string, $command_out, $command_log) = @$command;
	
	# check
	if (length($command_out) and length($command_log)) {
		# both 
		if (-e $command_out and -e $command_log) {
			print "=== Job: $command_string\n    previously finished, have $command_out and $command_log files\n" if $talk;
			return 1;
		}
		elsif (not -e $command_out and -e $command_log and 
			index($command_string, $opts{bamdedup}) == 0) 
		{
			# the deduplication command will not write out a bam file if the actual 
			# duplication rate is below the target rate
			# presume this is good!?
			print "=== Job: $command_string\n    presumed finished, have $command_log file only\n" if $talk;
			return 2;
		}
		elsif (-e $command_out and not -e $command_log) {
			# we have a output file but not a log file
			print "=== Job: $command_string\n    presumed finished, have $command_out file only, no log\n" if $talk;
			return 3;
		}
	}
	elsif (length($command_out)) {
		if (-e $command_out) {
			print "=== Job: $command_string\n    previously finished, have $command_out\n" if $talk;
			return 4;
		}
	}
	elsif (substr($command_string, 0, 2) eq 'rm') {
		# remove command doesn't leave an output (duh!) or log file
		# gotta check each one
		my $check = 0; 
		foreach my $item (split /\s+/, $command_string) {
			next if $item eq 'rm';
			$check++ if -e $item; # check true if file is present
		}
		if ($check == 0) {
			print "=== Job: $command_string\n    previously finished, target files missing\n" if $talk;
			return 5;
		}
	}
	
	# else presume command was not finished
	return 0; 
}

sub update_progress_file {
	my $key = shift;
	$progress{$key} = 1;
	return 1 if $opts{dryrun}; # just pretend
	my $fh = IO::File->new($progress_file, '>>') or 
		die "can't write to progress file! $!\n";
	$fh->print("$key\n");
	$fh->close;
	return 1;
}

sub run_input_peak_detection {
	return unless ($opts{blacklist} =~ /^(?:input|control)$/i);
	print "\n\n======= Generating exclusion list from reference control\n";
	if ($progress{control_peak}) {
		print "\nStep is completed\n";
		return;
	}
	
	# available reference bam files
	my @refbams;
	foreach my $Job (@Jobs) {
		if ( defined $Job->{control_bams} and scalar( @{$Job->{control_bams}} ) ) {
			push @refbams, @{$Job->{control_bams}};
		}
	}
	if (@refbams) {
		# we have at least one reference bam file to process
		# set the name of the exclusion list file
		$opts{blacklist} = File::Spec->catfile($opts{dir}, 
			sprintf("%s.control_peak", $opts{out}) );
	}
	else {
		# no reference bam files!
		print " No control reference bam files to process. Skipping\n";
		$opts{blacklist} = '';
		update_progress_file('control_peak');
		return;
	}
	
	# generate bam2wig command
	# very little filtering here - we basically want everything
	my $command = sprintf(
		"%s --out %s.bdg --nosecondary --noduplicate --nosupplementary --mean --bdg --bin %s --cpu %s ", 
		$opts{bam2wig}, 
		$opts{blacklist},
		$opts{chipbin},
		$opts{cpu} * $opts{job}, # give it everything we've got, single job
	);
	if ($opts{paired}) {
		$command .= sprintf("--span --pe --minsize %s --maxsize %s ", 
			$opts{minsize}, $opts{maxsize});
	}
	else {
		$command .= sprintf("--extend --extval %s ", $opts{fragsize});
	}
	if ($opts{chrskip}) {
		$command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
	}
	$command .= sprintf("--in %s ", join(',', @refbams));
	my $logfile = sprintf("%s.out.txt", $opts{blacklist});
	$command .= " 2>&1 > $logfile ";
	
	# add the mean bedgraph file
	$command .= sprintf(" && %s %s.bdg 2>> $logfile ", $opts{meanbdg}, $opts{blacklist});
	
	# add the q-value conversion
	$command .= sprintf(" && %s bdgcmp -t %s.bdg -c %s.global_mean.bdg -m qpois -o %s.qvalue.bdg 2>> $logfile ",
		$opts{macs}, $opts{blacklist}, $opts{blacklist}, $opts{blacklist});
	
	# add the peak call
	# we are using hard coded parameters for now, but I think these should be generic enough
	$command .= sprintf(" && %s bdgpeakcall -i %s.qvalue.bdg -c 3 -l 250 -g 500 --no-trackline -o %s.narrowPeak 2>> $logfile",
		$opts{macs}, $opts{blacklist}, $opts{blacklist});
	
	# add the peak conversion 
	$command .= sprintf(" && %s %s.narrowPeak 2>&1 >> $logfile", $opts{peak2bed}, 
		$opts{blacklist});
	
	# clean up
	$command .= sprintf(" && rm -f %s.bdg %s.global_mean.bdg %s.qvalue.bdg %s.narrowPeak %s.summit.bed",
		$opts{blacklist}, $opts{blacklist}, $opts{blacklist}, $opts{blacklist}, 
		$opts{blacklist});
	$opts{blacklist} .= '.bed'; # the actual output file
	
	# execute job
	execute_commands( [ [$command, $opts{blacklist}, $logfile], ] );
	update_progress_file('control_peak');
}

sub run_dedup {
	return unless ($opts{dedup});
	print "\n\n======= De-duplicating bam files\n";
	if ($progress{deduplication}) {
		print "\nStep is completed\n";
		return;
	}
	
	### Run de-duplication
	my @commands;
	my %name2done;
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_dedup_commands(\%name2done);
	}
	if (@commands) {
		execute_commands(\@commands);
		return if $opts{dryrun};
	}
	else {
		update_progress_file('deduplication');
		return;
	}
	
	### Collect deduplication statistics
	my @dedupstats;
	push @dedupstats, join("\t", qw(File TotalCount OpticalDuplicateCount DuplicateCount
		NonDuplicateCount DuplicationRate RetainedDuplicateCount));
	foreach my $c (@commands) {
		# initialize counts
		my $total   = 0; # total count
		my $optdup  = 0; # optical duplicate count
		my $nondup  = 0; # non-duplicate count
		my $dup     = 0; # duplicate count
		my $retdup  = 0; # retained duplicate count
		my $duprate = 0; # duplication rate
	
		# open log file and collect stats
		next if (not -e $c->[2]);
		my $fh = IO::File->new($c->[2], 'r');
		while (my $line = $fh->getline) {
			if ($line =~  /^  Total mapped:\s+(\d+)$/) {
				# bam_partial_dedup
				$total = $1;
			}
			elsif ($line =~ /^  Non-duplicate count:\s+(\d+)$/) {
				$nondup = $1;
			}
			elsif ($line =~ /^  Optical duplicate count:\s+(\d+)$/) {
				# optical bam_partial_dedup
				$optdup = $1;
			}
			elsif ($line =~ /^  Non-optical duplicate count:\s+(\d+)$/) {
				# non-optical bam_partial_dedup
				$dup = $1;
			}
			elsif ($line =~ /^  Non-optical duplication rate:\s+(\d\.\d+)$/) {
				# non-optical bam_partial_dedup
				$duprate = $1;
			}
			elsif ($line =~ /^  Retained non-optical duplicate count:\s+(\d+)\s*$/) {
				# bam_partial_dedup
				# oops, there may be a space at the end
				# this might not be found if no deduplication occurred
				$retdup = $1;
			}
		}
		$fh->close;
	
		# name of bam file, extracted from log file
		my (undef, undef, $name) = File::Spec->splitpath($c->[2]);
		$name =~ s/\.dedup\.out\.txt$//;
	
		# store in array
		push @dedupstats, join("\t", $name, $total, $optdup, $dup, $nondup, $duprate, 
			$retdup || $dup);
	}

	# print duplicate stats file
	my $dedupfile = File::Spec->catfile($opts{dir}, $opts{out} . '.dedup-stats.txt');
	my $fh = IO::File->new($dedupfile, 'w');
	foreach (@dedupstats) {
		$fh->print("$_\n");
	}
	$fh->close;
	print "\n Wrote deduplication report $dedupfile\n";
	
	update_progress_file('deduplication');
}

sub run_bam_check {
	print "\n\n======= Checking bam files\n";
	foreach my $Job (@Jobs) {
		$Job->find_dedup_bams;
	}
}

sub run_mappable_space_report {
	print "\n\n======= Determining mappable genome space\n";
	
	# check the user supplied value
	if ($opts{genome} and not $opts{dryrun} and not $progress{mappable_size}) {
		
		# get the full genome size from the chromosome file
		my $fh = IO::File->new($chromofile, 'r');
		my $genome_size = 0;
		while (my $line = $fh->getline) {
			if ($line =~ /\s+(\d+)$/) {
				$genome_size += $1;
			}
		}
		$fh->close;
		
		# double-check the user value
		my $ratio = $opts{genome} / $genome_size;
		if ($ratio > 1) {
			printf "\nUser supplied genome size (%d) is larger than actual genome size (%d)!!!\nDetermining actual empirical mappable size\n",
				$opts{genome}, $genome_size;
			$opts{genome} = 0;
		}
		elsif ($ratio < 0.6) {
			printf "\nUser supplied genome size (%d) is considerably smaller than actual genome size (%d)!\nDetermining actual empirical mappable size\n",
				$opts{genome}, $genome_size;
			$opts{genome} = 0;
		}
		else {
			# ratio is somewhere between 50-100% of actual genome size, so assume ok
			printf "\nUsing user-specified size of %d\n", $opts{genome};
			return;
		}
	}
	elsif ($opts{genome} and $opts{dryrun}) {
		printf "\nPretending that user supplied genome size is ok\n";
		return;
	}
	
	# the output logfile 
	my $logfile = File::Spec->catfile($opts{dir}, 
		sprintf("%s.mappable.out.txt", $opts{out}) );
	
	# check if command is finished, otherwise run it
	if ($progress{mappable_size}) {
		print "\nStep is completed\n";
		# we will read the value below
	}
	else {
		# collect all available bam files
		my @bamlist;
		foreach my $Job (@Jobs) {
			push @bamlist, @{$Job->{control_use_bams}};
			push @bamlist, @{$Job->{chip_use_bams}}
		}
		my $command = sprintf("%s --cpu %d ", $opts{reportmap}, $opts{cpu} * $opts{job});
			# give the cpu everything we've got, there's only one job
		if ($opts{chrskip}) {
			$command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		$command .= join(" ", @bamlist);
		$command .= " 2>&1 > $logfile";
		
		# execute
		# the log file is the output
		execute_commands( [ [$command, $logfile, $logfile] ] );
	}
	
	# Collect results from output file
	# do this regardless whether this was finished previously, since we have to 
	# extract the value into memory anyway
	if (-e $logfile) {
		my $fh = IO::File->new($logfile, 'r') or 
			die " unable to open mappable report file '$logfile'!";
		while (my $line = $fh->getline) {
			# we're going to use the all mappable space number
			if ($line =~ /All mappable space: ([\d\.]+) Mb/) {
				$opts{genome} = $1 * 1000000;
				last;
			}
		}
		$fh->close;
		if ($opts{genome}) {
			printf "\n Genome mappable space calculated to be %d bp\n", $opts{genome};
		}
		else {
			die "\n Unable to extract genome mappable space from report '$logfile'!\n";
		}
	}
	elsif ($opts{dryrun}) {
		print "\n Artificially setting mappable genome size to 100000000 (100 Mb) for dry run purposes\n";
	}
	else {
		die "\n Genome mappable report log file is not present! Unable to continue!\n";
	}
	
	update_progress_file('mappable_size'); # might be redundant but that's ok
}


sub run_bam_fragment_conversion {
	my @commands;
	my %name2done;
	print "\n\n======= Generating fragment coverage files\n";
	
	
	# Generate commands for each job
	# we need the command regardless of whether it needs to be run or not
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_bam2wig_frag_commands(\%name2done);
	}
	
	# Execute as necessary
	if ($progress{fragment}) {
		print "\nStep is completed\n";
	}
	elsif (@commands) {
		# run programs
		execute_commands(\@commands);
		update_progress_file('fragment');
	}
	
	# skip counting results if dryrun
	if ($opts{dryrun}) {
		# artificially set target depth
		print "\n Artificially setting target depth to 25 Million for dry run purposes\n";
		$opts{targetdep} = 25;
		return;
	}
	
	# count results
	my %bam2count;
	foreach my $com (@commands) {
		my $log = $com->[2];
		my $fh = IO::File->new($log, 'r') or  # this should open!!!!
			die "something is wrong! Job completed but unable to open $log!? $!";
		
		my @files; # there may be one or more files processed here
		while (my $line = $fh->getline) {
			chomp $line;
			if ($line =~ /^ Processing files (.+)\.\.\.$/) {
				# the names of files
				@files = split(/, /, $1);
			}
			elsif ($line =~ /^ Normalizing depth based on ([\d,]+) total counted (?:alignments|fragments)$/) {
				# only one file was processed
				# need to grab the name from the list
				my $count = $1;
				unless (exists $bam2count{$files[0]}) {
					$count =~ s/,//g; # remove commas
					$bam2count{$files[0]} = sprintf "%.1f", $count / 1_000_000;
				}
			}
			elsif ($line =~ /^  (.+) had ([\d,]+) total counted (?:alignments|fragments)$/) {
				# multiple files were processed
				my $file = $1;
				my $count = $2;
				unless (exists $bam2count{$file}) {
					$count =~ s/,//g; # remove commas
					$bam2count{$file} = sprintf "%.1f", $count / 1_000_000;
				}
			}
		}
		$fh->close;
	}
	
	# print report
	print "\n Total fragments accepted for analysis\n";
	foreach my $f (sort {$a cmp $b} keys %bam2count) {
		printf "  %6sM  $f\n", $bam2count{$f};
	}
	
	# Calculate minimum target depth to use
	my $targetdep = int( min(values %bam2count) ) || 1; # just in case!
	if (defined $opts{targetdep}) {
		printf "\n WARNING!!! Calculated target sequence depth of %d Million is overridden by manually set value %d\n",
			$targetdep, $opts{targetdep};
	}
	else {
		$opts{targetdep} = $targetdep;
		printf "\n Setting target sequence depth to %d Million\n", $opts{targetdep};
	}
}

sub run_bam_count_conversion {
	my @commands;
	my %name2done;
	print "\n\n======= Generating fragment count files\n";
	if ($progress{count}) {
		print "\nStep is completed\n";
		return;
	}
	
	# generate commands and run
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_bam2wig_count_commands(\%name2done);
	}
	if (@commands) {
		# run programs
		execute_commands(\@commands);
	}
	
	update_progress_file('count');
}

sub run_lambda_control {
	my @commands;
	my %name2done;
	print "\n\n======= Generating control files\n";
	if ($progress{lambda}) {
		print "\nStep is completed\n";
		return;
	}
	foreach my $Job (@Jobs) {
		# this handles either lambda_control or global mean files
		push @commands, $Job->generate_lambda_control_commands(\%name2done);
	}
	if (@commands) {
		execute_commands(\@commands);
	}
	update_progress_file('lambda');
}

sub run_bw_conversion {
	my @commands;
	my %name2done;
	print "\n\n======= Converting Fragment bigWig files to bedGraph\n";
	if ($progress{bw2bdg}) {
		print "\nStep is completed\n";
		return;
	}
	foreach my $Job (@Jobs) {
		push @commands, $Job->convert_bw_to_bdg(\%name2done);
	}
	if (@commands) {
		execute_commands(\@commands);
	}
	update_progress_file('bw2bdg');
}

sub run_bdgcmp {
	my @commands;
	print "\n\n======= Generate enrichment files\n";
	if ($progress{bdgcmp}) {
		print "\nStep is completed\n";
		return;
	}
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_enrichment_commands();
	}
	execute_commands(\@commands);
	update_progress_file('bdgcmp');
}

sub run_call_peaks {
	my @commands;
	print "\n\n======= Call peaks\n";
	if ($progress{callpeak}) {
		print "\nStep is completed\n";
		return;
	}
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_peakcall_commands;
	}
	execute_commands(\@commands);
	update_progress_file('callpeak');
}

sub run_clean_peaks {
	my @commands;
	print "\n\n======= Cleaning peak files\n";
	if ($progress{cleanpeak}) {
		print "\nStep is completed\n";
		return;
	}
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_cleanpeak_commands;
	}
	execute_commands(\@commands);
	update_progress_file('cleanpeak');
}

sub run_bdg_conversion {
	my @commands;
	my %name2done;
	print "\n\n======= Converting bedGraph files\n";
	if ($progress{bdg2bw}) {
		print "\nStep is completed\n";
		return;
	}
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_bdg2bw_commands($chromofile, \%name2done);
	}
	execute_commands(\@commands);
	update_progress_file('bdg2bw');
}

sub run_peak_merge {
	return if scalar(@Jobs) == 1; # no sense merging one job!
	print "\n\n======= Merging called Peak files\n";
	if ($progress{peakmerge}) {
		print "\nStep is completed\n";
		return;
	}
	die "no bedtools application in path!\n" unless $opts{bedtools} =~ /\w+/;
	die "no intersect_peaks.pl application in path!\n" unless $opts{intersect} =~ /\w+/;
	
	my @commands;
	
	# narrowPeaks
	my $merge_file = File::Spec->catfile($opts{dir}, $opts{out});
	my $command = sprintf("%s --name %s_merge --bed %s --out %s --genome %s ", 
		$opts{intersect}, $opts{out}, $opts{bedtools}, $merge_file, $chromofile);
	my $command_check = length($command);
	my $count_check = 0;
	foreach my $Job (@Jobs) {
		if (
			$Job->{clean_peak} and 
			( (-e $Job->{clean_peak} and -s _ > 0) or $opts{dryrun} )
		) {
			# only add the file if it exists and non-zero in length
			# or just fake it if we're running a dry run
			$command .= sprintf("%s ", $Job->{clean_peak});
			$count_check++;
		}
	}
	if (length($command) > $command_check and $count_check > 1) {
		my $log = $merge_file . '.merge.out.txt';
		$command .= sprintf("2>&1 > %s ", $log);
		push @commands, [$command, $merge_file . '.bed', $log]; 
			# this will have multiple outputs, but one is just a .bed file
	}
	else {
		unless ($opts{dryrun}) {
			print "One or fewer narrow peak files found, nothing to merge!\n";
		}
	}
	
	# broadPeaks
	if ($opts{broad}) {
		my $merge2_file = File::Spec->catfile($opts{dir}, $opts{out} . "_broad");
		my $command2 = sprintf("%s --name %s_gapmerge --bed %s --out %s --genome %s ", 
			$opts{intersect}, $opts{out}, $opts{bedtools}, $merge2_file, $chromofile);
		my $command2_check = length($command2);
		my $count2_check = 0;
		
		foreach my $Job (@Jobs) {
			if (
				$Job->{clean_gappeak} and 
				( (-e $Job->{clean_gappeak} and -s _  > 0) or $opts{dryrun} )
			) {
				# only add the file if it exists and non-zero in length
				# or just fake it if we're running a dry run
				$command2 .= sprintf("%s ", $Job->{clean_gappeak});
				$count2_check++;
			}
		}
		if (length($command2) > $command2_check and $count2_check > 1) {
			my $log2 = $merge2_file . '.merge.out.txt';
			$command2 .= sprintf("2>&1 > %s ", $log2);
			push @commands, [$command2, $merge2_file . '.bed', $log2];
				# this will have multiple outputs, but one is just a .bed file
		}
		else {
			unless ($opts{dryrun}) {
				print "One or fewer gapped peak files found, nothing to merge!\n";
			}
		}
	}
	
	if (@commands) {
		execute_commands(\@commands);
		update_progress_file('peakmerge');
	}
}

sub run_rescore {
	print "\n\n======= Re-scoring all peaks\n";
	if ($progress{rescore}) {
		print "\nStep is completed\n";
		return;
	}
	die "no get_datasets.pl script in path!\n" unless $opts{getdata} =~ /\w+/;
	die "no get_relative_data.pl script in path!\n" unless $opts{getrel} =~ /\w+/;
	
	# prepare filenames
	my $input;
	if (scalar(@Jobs) > 1) {
		$input = File::Spec->catfile($opts{dir}, $opts{out} . '.bed');
		if (not $opts{dryrun} and not -e $input) {
			print "No merged peak file '$input'!\n";
			return;
		}
	}
	else {
		$input = $Jobs[0]->{clean_peak};
		if (not $opts{dryrun} and not -e $input) {
			print "No peak file '$input'!\n";
			return;
		}
	}
	my $output1 = File::Spec->catfile($opts{dir}, $opts{out} . '_meanQvalue.txt');
	my $output2 = File::Spec->catfile($opts{dir}, $opts{out} . '_meanLog2FE.txt');
	my $output3 = File::Spec->catfile($opts{dir}, $opts{out} . '_counts.txt');
	my $output4 = File::Spec->catfile($opts{dir}, $opts{out} . '_profile_fragment.txt');
	my $output5 = File::Spec->catfile($opts{dir}, $opts{out} . '_profile_log2FE.txt');
	my $output6 = File::Spec->catfile($opts{dir}, $opts{out} . '_genome_counts.txt.gz');
	
	# start list of conditions
	my @conditions = ("Replicate\tDataset\n");
	
	# generate four get_dataset and two get_relative commands
	# go ahead and make the fourth genome-wide command, even though we may not use it
	my $command1 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
		$opts{getdata}, $opts{cpu}, $input, $output1);
	my $command2 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
		$opts{getdata}, $opts{cpu}, $input, $output2);
	my $command3 = sprintf("%s --method sum --cpu %s --in %s --out %s --format 0 ",
		$opts{getdata}, $opts{cpu}, $input, $output3);
	my $command4 = sprintf("%s --method mean --cpu %s --in %s --out %s --win %s --num 25 --pos m --long --format 3 --groups --sum ",
		$opts{getrel}, $opts{cpu}, $input, $output4, $opts{binsize});
	my $command5 = sprintf("%s --method mean --cpu %s --in %s --out %s --win %s --num 25 --pos m --long --format 3 --groups --sum ",
		$opts{getrel}, $opts{cpu}, $input, $output5, $opts{binsize});
	my $command6 = sprintf("%s --method sum --cpu %s --feature genome --win %d --discard %s --out %s --format 0 ",
		$opts{getdata}, $opts{cpu}, $opts{window}, $opts{discard}, $output6);
	
	# add dataset files
	my %name2done;
	foreach my $Job (@Jobs) {
		if ($Job->{qvalue_bw}) {
			$command1 .= sprintf("--data %s ", $Job->{qvalue_bw});
		}
		if ($Job->{logfe_bw}) {
			$command2 .= sprintf("--data %s ", $Job->{logfe_bw});
			$command5 .= sprintf("--data %s ", $Job->{logfe_bw});
		}
		if ($Job->{chip_bw}) {
			$command4 .= sprintf("--data %s ", $Job->{chip_bw});
		}
		foreach my $b ( @{ $Job->{chip_count_bw} } ) {
			$command3 .=  "--data $b ";
			$command6 .=  "--data $b ";
			my $name = simplify_dataset_name($b);
			push @conditions, sprintf("%s\t%s\n", $name, $Job->{name})
		}
		foreach my $b ( @{ $Job->{control_count_bw} } ) {
			next if exists $name2done{$b};
			$command3 .=  "--data $b ";
			$command6 .=  "--data $b ";
			my $name = simplify_dataset_name($b);
			push @conditions, "$name\tInput\n";
			$name2done{$b} = 1; # remember it's done
		}
	}
	
	
	# add log outputs to commands
	my @commands;
	my $log = $output1;
	$log =~ s/txt$/out.txt/;
	$command1 .= " 2>&1 > $log";
	push @commands, [$command1, $output1, $log];
	
	$log = $output2;
	$log =~ s/txt$/out.txt/;
	$command2 .= " 2>&1 > $log";
	push @commands, [$command2, $output2, $log];
	
	$log = $output3;
	$log =~ s/txt$/out.txt/;
	$command3 .= " 2>&1 > $log";
	push @commands, [$command3, $output3, $log];
	
	$log = $output4;
	$log =~ s/txt$/out.txt/;
	$command4 .= " 2>&1 > $log";
	push @commands, [$command4, $output4, $log];
	
	$log = $output5;
	$log =~ s/txt$/out.txt/;
	$command5 .= " 2>&1 > $log";
	push @commands, [$command5, $output5, $log];
	
	if ($opts{window}) {
		# user has given an actual genome window size, so we'll run this command
		$log = $output6;
		$log =~ s/txt\.gz$/out.txt/;
		$command6 .= " 2>&1 > $log";
		push @commands, [$command6, $output6, $log];
	}
	
	
	# broad peak rescore
	my ($output7, $output8, $output9);
	if ($opts{broad}) {
		# generate broad file names
		my $input2;
		if (scalar(@Jobs) > 1) {
			$input2 = File::Spec->catfile($opts{dir}, $opts{out} . '_broad.bed');
			if (not $opts{dryrun} and not -e $input2) {
				die "unable to find merged gapped Peak bed file '$input2'!\n";
			}
		}
		else {
			$input2 = $Jobs[0]->{clean_gappeak};
			if (not $opts{dryrun} and not -e $input2) {
				die "unable to find gapped Peak bed file '$input2'!\n";
			}
		}
		$output7 = File::Spec->catfile($opts{dir}, $opts{out} . '_broad_meanQvalue.txt');
		$output8 = File::Spec->catfile($opts{dir}, $opts{out} . '_broad_meanLog2FE.txt');
		$output9 = File::Spec->catfile($opts{dir}, $opts{out} . '_broad_counts.txt');
		
		# generate three get_dataset commands
		# we won't run get_relative data for gapped peaks
		my $command7 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
			$opts{getdata}, $opts{cpu}, $input2, $output7);
		my $command8 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
			$opts{getdata}, $opts{cpu}, $input2, $output8);
		my $command9 = sprintf("%s --method sum --cpu %s --in %s --out %s --format 0 ",
			$opts{getdata}, $opts{cpu}, $input2, $output9);
		%name2done = ();
		foreach my $Job (@Jobs) {
			if ($Job->{qvalue_bw}) {
				$command7 .= sprintf("--data %s ", $Job->{qvalue_bw});
			}
			if ($Job->{logfe_bw}) {
				$command8 .= sprintf("--data %s ", $Job->{logfe_bw});
			}
			foreach my $b ( @{ $Job->{chip_count_bw} } ) {
				$command9 .=  "--data $b ";
			}
			foreach my $b ( @{ $Job->{control_count_bw} } ) {
				next if exists $name2done{$b};
				$command9 .=  "--data $b ";
				$name2done{$b} = 1; # remember it's done
			}
		}
	
		# add log outputs to commands
		my @commands;
		my $log = $output7;
		$log =~ s/txt$/out.txt/;
		$command7 .= " 2>&1 > $log";
		push @commands, [$command7, $output7, $log];
		
		$log = $output8;
		$log =~ s/txt$/out.txt/;
		$command8 .= " 2>&1 > $log";
		push @commands, [$command8, $output8, $log];
		
		$log = $output9;
		$log =~ s/txt$/out.txt/;
		$command9 .= " 2>&1 > $log";
		push @commands, [$command9, $output9, $log];
	}
	
	# write conditions file
	my $samplefile = File::Spec->catfile($opts{dir}, $opts{out} . '_samples.txt');
	unless ($opts{dryrun}) {
		my $fh = IO::File->new($samplefile, "w");
		foreach (@conditions) {
			$fh->print($_);
		}
		$fh->close;
	}
	
	
	
	### Execute commands
	execute_commands(\@commands);
	
	### Replicate Merge
	# do this here so that we still know the output commands
	# must be done after the data collection anyway, so can't be merged above
	if ($opts{repmean}) {
		# we will generate count means of the replicates
		return unless $opts{combrep} =~ /\w+/;
		print "\n\n======= Generating sample replicate means\n";
		my @commands2;
		
		# narrowPeak counts
		my $output3m = File::Spec->catfile($opts{dir}, $opts{out} . '_meanCounts.txt');
		my $command10 = sprintf("%s --in %s --out %s --sample %s --method mean --format 0 ", 
			$opts{combrep}, $output3, $output3m, $samplefile);
		$log = $output3m;
		$log =~ s/txt/out.txt/;
		$command10 .= " 2>&1 > $log";
		push @commands2, [$command10, $output3m, $log];
		
		# genome counts
		if ($opts{window}) {
			my $output6m = File::Spec->catfile($opts{dir}, $opts{out} . 
				'_genome_meanCounts.txt.gz');
			my $command11 = sprintf("%s --in %s --out %s --sample %s --method mean --format 0 ", 
				$opts{combrep}, $output6, $output6m, $samplefile);
			$log = $output6m;
			$log =~ s/txt\.gz/out.txt/;
			$command11 .= " 2>&1 > $log";
			push @commands2, [$command11, $output6m, $log];
		}
		
		# broadPeak Counts
		if ($opts{broad}) {
			my $output9m = File::Spec->catfile($opts{dir}, $opts{out} . 
				'_broad_meanCounts.txt');
			my $command12 = sprintf("%s --in %s --out %s --sample %s --method mean --format 0 ", 
				$opts{combrep}, $output9, $output9m, $samplefile);
			$log = $output9m;
			$log =~ s/txt/out.txt/;
			$command12 .= " 2>&1 > $log";
			push @commands2, [$command12, $output9m, $log];
		}
		
		# execute
		execute_commands(\@commands2);
	}
	
	update_progress_file('rescore');
}

sub run_efficiency {
	print "\n\n======= Scoring peak calls for ChIP efficiency\n";
	if ($progress{efficiency}) {
		print "\nStep is completed\n";
		return;
	}
	die "no get_chip_efficiency.pl script in path!\n" unless $opts{geteff} =~ /\w+/;
	
	### Sample conditions file
	my $samplefile = File::Spec->catfile($opts{dir}, $opts{out} . '_samples.txt');
	if (not $opts{dryrun} and not -e $samplefile) {
		print "No sample file '$samplefile'!\n";
		return;
	}
	
	### ChIP efficiency
	my @commands;
	
	# universal control counts
	# I have to search for these explicitly, since they're not associated with ChIP job
	my @universal_counts;
	foreach my $Job (@Jobs) {
		next unless (
			not defined $Job->{clean_peak} and 
			scalar @{ $Job->{control_count_bw} } > 0
		);
		push @universal_counts, @{ $Job->{control_count_bw} };
	}
	
	# collect count files for each job
	foreach my $Job (@Jobs) {
		next unless defined $Job->{clean_peak};
		my $output = File::Spec->catfile($opts{dir}, $Job->{name} . '.efficiency.txt');
		my $command = sprintf("%s --in %s --group %s --out %s --cpu %d ", 
			$opts{geteff}, $Job->{clean_peak}, $samplefile, $output, $opts{cpu});
		# add count files, we should have at least one chip and one control
		foreach my $b ( @{ $Job->{chip_count_bw} } ) {
			$command .= "$b ";
		}
		foreach my $b ( @{ $Job->{control_count_bw} } ) {
			$command .= "$b ";
		}
		foreach my $b (@universal_counts) {
			$command .= "$b ";
		}
		my $log = $output;
		$log =~ s/txt/out.txt/;
		$command .= sprintf("2>&1 > $log");
		push @commands, [$command, $output, $log];
	}
	
	# execute the efficiency commands
	execute_commands(\@commands);
	
	# proceed no further if dry run
	return if $opts{dryrun};
	
	# merge the efficiency outputs into one
	if (scalar @commands > 1) {
		my @combined_eff_data;
		my @combined_eff_meta;
		my $eff_header;
		foreach my $c (@commands) {
			my $fh = IO::File->new($c->[1], 'r');
			while (my $line = $fh->getline) {
				if (substr($line, 0, 1) eq '#') {
					push @combined_eff_meta, $line;
				}
				elsif ($line =~ /^Replicate/) {
					$eff_header = $line unless $eff_header;
				}
				else {
					push @combined_eff_data, $line;
				}
			}
			$fh->close;
			unlink $c->[1];
		}
		
		# write merged file
		my $combined_eff_out = File::Spec->catfile($opts{dir}, $opts{out} . '.chip_efficiency.txt');
		my $fh = IO::File->new($combined_eff_out, 'w');
		foreach (@combined_eff_meta) {
			$fh->print($_);
		}
		$fh->print($eff_header);
		foreach (@combined_eff_data) {
			$fh->print($_);
		}
		$fh->close;
		print "\nWrote combined ChIP efficiency file $combined_eff_out\n";
	}
	elsif (scalar @commands == 1) {
		# just one, so rename it for consistency sake
		my $eff_out = File::Spec->catfile($opts{dir}, $opts{out} . '.chip_efficiency.txt');
		move($commands[0]->[1], $eff_out);
	}
	
	update_progress_file('efficiency');
}

sub run_plot_peaks {
	return unless ($opts{plot});
	print "\n\n======= Plotting Peak figures\n";
	if ($opts{rscript} !~ /\w+/ or $opts{plotpeak} !~ /\w+/) {
		print "Rscript or plot_peak_figures.R script not defined!\n";
		return;
	}
	my $outbase = File::Spec->catfile($opts{dir}, $opts{out});
	if (not $opts{dryrun} and not -e "$outbase.bed") {
		print "No output files for plotting\n";
		return;
	}
	my $command = sprintf("%s %s --input %s ", $opts{rscript}, $opts{plotpeak}, $outbase);
	my $log = $outbase . '_plot_figures.out.txt';
	$command .= " 2>&1 > $log";
	
	# there are multiple output files from this script
	# only using one as an example
	my $example = File::Spec->catfile($opts{dir}, $opts{out} . '_PCA.png');
	execute_commands( [ [$command, $example, $log] ] );
}


sub run_cleanup {
	return if $opts{dryrun};
	print "\n\n======= Combining log files\n";
	
	
	# combine output logs
	my @combined_output;
	my @logs = glob(File::Spec->catfile($opts{dir}, '*.out.txt'));
	foreach my $log (@logs) {
		push @combined_output, "=== Log file: $log\n";
		if (-z $log) {
			# an empty file
			push @combined_output, "\n";
			unlink $log;
		}
		else {
			# push log contents to combined output
			my $fh = IO::File->new($log, 'r') or next;
			push @combined_output, <$fh>;
			push @combined_output, "\n";
			$fh->close;
		}
		unlink $log if -e $log;
	}
	my $file = File::Spec->catfile($opts{dir}, $opts{out} . "_job_output_logs.txt");
	my $fh = IO::File->new($file, "w");
	foreach (@combined_output) {
		$fh->print($_);
	}
	$fh->close;
	print "\nCombined all job output log files into '$file'\n";
	
	# remove files no longer need
	print "\n\n======= Deleting temporary files\n";
	unless ($opts{savebam}) {
		foreach my $Job (@Jobs) {
			foreach my $b ( @{ $Job->{chip_dedup_bams} } ) {
				unlink $b if -e $b;
				$b .= ".bai";
				unlink $b if -e $b;
			}
			foreach my $b ( @{ $Job->{control_dedup_bams} } ) {
				unlink $b if -e $b;
				$b .= ".bai";
				unlink $b if -e $b;
			}
		}
	}
	unlink $chromofile;
	unlink $progress_file;
}


sub run_organize {
	return unless ($opts{organize});
	return if $opts{dryrun};
	print "\n\n======= Moving files into subdirectories\n";
	
	# directories
	my $fragdir  = File::Spec->catfile($opts{dir}, 'Fragment');
	my $log2dir  = File::Spec->catfile($opts{dir}, 'Log2FE');
	my $countdir = File::Spec->catfile($opts{dir}, 'Count');
	my $qdir     = File::Spec->catfile($opts{dir}, 'QValue');
	my $peakdir  = File::Spec->catfile($opts{dir}, 'Peaks');
	my $sumitdir = File::Spec->catfile($opts{dir}, 'PeakSummits');
	my $analdir  = File::Spec->catfile($opts{dir}, 'Analysis');
	foreach ($fragdir, $log2dir, $countdir, $qdir, $peakdir, $sumitdir, $analdir) {
		make_path($_);
	}
	
	# we're globbing all the files, hope this isn't a problem in case user has 
	# similarly named files in existing directory.
	# Otherwise there's an awful lot of conditional checks for files in every single 
	# ChIPJob object plus general run files....
	
	# log2FE files
	foreach (glob(File::Spec->catfile($opts{dir}, '*.log2FE.bw')) ) {
		move($_, $log2dir);
	}
	
	# fragment files
	foreach (glob(File::Spec->catfile($opts{dir}, '*.fragment.bw')) ) {
		move($_, $fragdir);
	}
	foreach (glob(File::Spec->catfile($opts{dir}, '*.lambda_control.bw')) ) {
		move($_, $fragdir);
	}
	foreach (glob(File::Spec->catfile($opts{dir}, '*.control_fragment.bw')) ) {
		move($_, $fragdir);
	}
	foreach (glob(File::Spec->catfile($opts{dir}, '*.fragment.global_mean.bw')) ) {
		move($_, $fragdir);
	}
	
	# qvalue files
	foreach (glob(File::Spec->catfile($opts{dir}, '*.qvalue.bw')) ) {
		move($_, $qdir);
	}
	
	# count files
	foreach (glob(File::Spec->catfile($opts{dir}, '*.count.bw')) ) {
		move($_, $countdir);
	}
	
	# merged peak
	foreach (glob(File::Spec->catfile($opts{dir}, '*.summit.bed')) ) {
		move($_, $sumitdir);
	}
	foreach (glob(File::Spec->catfile($opts{dir}, '*.bed')) ) {
		move($_, $peakdir);
	}
	
	# text files
	foreach (glob(File::Spec->catfile($opts{dir}, '*.txt*')) ) {
		next if $_ =~ /job_output_logs\.txt$/;
		move($_, $analdir);
	}
	
	# image files
	if ($opts{plot}) {
		my $imagedir = File::Spec->catfile($opts{dir}, 'Plots');
		make_path($imagedir);
		foreach (glob(File::Spec->catfile($opts{dir}, '*.png')) ) {
			move($_, $imagedir);
		}
	}
	
	# dedup bam files
	if ($opts{savebam} and $opts{dedup}) {
		my $bamdir = File::Spec->catfile($opts{dir}, 'DeDupBam');
		make_path($bamdir);
		foreach (glob(File::Spec->catfile($opts{dir}, '*.dedup.bam*')) ) {
			move($_, $bamdir);
		}
	}
	
	# saved bedGraph files
	if ($opts{savebdg}) {
		my $bdgdir = File::Spec->catfile($opts{dir}, 'BedGraph');
		make_path($bdgdir);
		foreach (glob(File::Spec->catfile($opts{dir}, '*.bdg')) ) {
			move($_, $bdgdir);
		}
	}
}


############### ChIPJob Package ########################################################

package ChIPjob;
use Carp;
use Data::Dumper;

sub new {
	# pass name, comma-list of ChIP files, and comma-list of control files
	my ($class, $name, $chip, $control, $chip_scale, $control_scale, $chrnorm) = @_;
	my $namepath = File::Spec->catfile($opts{dir}, $name);
	my $self = bless {
	    name                => $name, # name for this ChIP job
	    chip_bams           => undef, # array of initial chip bam file names
	    control_bams        => undef, # array of initial conrol bam file names
	    chip_dedup_bams     => [], # array of deduplicated bam file names
	    chip_use_bams       => [], # array of final bam file names to use
	    control_dedup_bams  => [], # array of deduplicated control bam files
	    control_use_bams    => [], # array of final control bam files to use
	    chip_scale          => undef, # array of scaling factors for chip signal
	    control_scale       => undef, # array of scaling factors for control signal
	    chrnorm             => $chrnorm, # chromosome-scaling factor
	    chip_bw             => undef, # ChIP signal bigWig file name
	    chip_bdg            => undef, # ChIP signal bedGraph file name
	    chip_count_bw       => [], # array of ChIP count bigWig file names
	    control_count_bw    => [], # array of countrol count bigWig file names
	    d_control_bdg       => undef, # d control pileup bedGraph
	    d_control_bw        => undef, # d control pileup bigWig file name
	    s_control_bdg       => undef, # small lambda control bedGraph file name
	    l_control_bdg       => undef, # large lambda control bedGraph file name
	    sld_control_file    => undef, # combined d, small, large lambda file
	    lambda_bdg          => undef, # final lambda control bedGraph file name
	    lambda_bw           => undef, # final lambda control bigWig file name
	    fe_bdg              => undef, # fold enrichment bedGraph file name
	    logfe_bw            => undef, # log2 fold enrichment bigWig file name
	    qvalue_bdg          => undef, # q-value bedGraph file name
	    qvalue_bw           => undef, # q-value bigWig file name
	    peak                => undef, # narrow peak file name
	    gappeak             => undef, # gapped peak file name
	    clean_peak          => undef, # cleaned peak file name
	    peak_summit         => undef, # cleaned peak summit file name
	    clean_gappeak       => undef, # cleaned gapped peak file name
	}, $class;
	
	
	## check the ChIP files
	if ($chip =~ /\.(?:bw|bigwig)$/i) {
		$self->{chip_bw} = $chip;
		$self->crash("only one ChIP bigWig file is allowed per experiment!\n") if 
			$chip =~ /,/;
		$self->{chip_bdg} = "$namepath.fragment.bdg";
		$self->{qvalue_bdg} = "$namepath.qvalue.bdg";
		$self->{qvalue_bw} = "$namepath.qvalue.bw";
		$self->{fe_bdg} = "$namepath.FE.bdg";
		$self->{logfe_bw} = "$namepath.log2FE.bw";
		$self->{peak} = "$namepath.narrowPeak";
		$self->{gappeak} = "$namepath.gappedPeak";
		$self->{clean_peak} = "$namepath.bed";
		$self->{peak_summit} = $namepath . '.summit.bed';
		$self->{clean_gappeak} = "$namepath.gapped.bed";
	}
	elsif ($chip =~ /\.bam$/i) {
		my @bams = split(',', $chip);
		$self->{chip_bams} = \@bams;
		$self->{chip_bdg} = "$namepath.fragment.bdg";
		$self->{chip_bw} = "$namepath.fragment.bw";
		$self->{qvalue_bdg} = "$namepath.qvalue.bdg";
		$self->{qvalue_bw} = "$namepath.qvalue.bw";
		$self->{fe_bdg} = "$namepath.FE.bdg";
		$self->{logfe_bw} = "$namepath.log2FE.bw";
		$self->{peak} = "$namepath.narrowPeak";
		$self->{gappeak} = "$namepath.gappedPeak";
		$self->{clean_peak} = "$namepath.bed";
		$self->{clean_gappeak} = "$namepath.gapped.bed";
		if ($chip_scale) {
			$self->{chip_scale} = [split(',', $chip_scale)];
			$self->crash("unequal scale factors and bam files!\n") if 
				scalar(@{$self->{chip_bams}}) != scalar(@bams);
		}
		# generate count bw and dedup bam file names
		foreach my $bam (@bams) {
			my (undef, undef, $fname) = File::Spec->splitpath($bam);
			$fname =~ s/\.bam$//i; # strip extension
			my $base = File::Spec->catfile($opts{dir}, $fname);
			# count bigWig
			push @{ $self->{chip_count_bw} }, "$base.count.bw";
			# dedup bam
			push @{ $self->{chip_dedup_bams} }, "$base.dedup.bam";
		}
	}
	else {
		# must be a control
		$self->{chip_bams} = [];
	}
	
	
	## check the control files
	if ($control =~ /^Custom\-Universal\-(.+)$/) {
		# this is just a ChIP Job using a common universal reference file
		# the actual control name is embedded in the name, so extract it
		# we only need to know the bedgraph file name here
		$self->{lambda_bdg} = File::Spec->catfile($opts{dir}, $1 . '.bdg'); 
	}
	elsif ($control =~ /\.(?:bw|bigwig)$/i) {
		# pre-generated reference bigWig file
		$self->crash("only one control bigWig file is allowed per experiment!\n") if
			$chip =~ /,/;
		$self->{lambda_bw} = $control;
		$self->{lambda_bdg} = $namepath . '.bdg'; 
	}
	elsif ($control =~ /\.bam$/i) {
		# we have control bams to process
		my @bams = split(',', $control);
		$self->{control_bams} = \@bams;
		if ($opts{lambda}) {
			my $control_base = $namepath;
			if ($namepath =~ /\.lambda_control/) {
				# the lambda_control extension is already appended to the name 
				# as is the case with universal controls
				$self->{lambda_bdg} = $namepath . '.bdg';
				$self->{lambda_bw} = $namepath . '.bw';
				# strip the lambda_control bit for the other intermediate files
				$control_base =~ s/\.lambda_control//;
			}
			else {
				# append the lambda_control extension
				$self->{lambda_bdg} = $namepath . '.lambda_control.bdg';
				$self->{lambda_bw} = $namepath . '.lambda_control.bw';
			}
			$self->{d_control_bdg} = "$control_base.dlocal.bdg";
			$self->{d_control_bw}  = "$control_base.control_fragment.bw";
			$self->{s_control_bdg} = "$control_base.slocal.bdg" if $opts{slocal};
			$self->{l_control_bdg} = "$control_base.llocal.bdg" if $opts{llocal};
			$self->{sld_control_file} = "$control_base.sldlocal.txt";
		}
		else {
			# even with no lambda, we recycle the hash entries for simplicity
			if ($namepath =~ /\.control_fragment/) {
				# the fragment control extension is already appended to the name
				# as is the case with universal controls
				$self->{lambda_bdg} = $namepath . '.bdg';
				$self->{lambda_bw} = $namepath . '.bw';
			}
			else {
				# append the fragment control extension
				$self->{lambda_bdg} = $namepath . '.control_fragment.bdg';
				$self->{lambda_bw} = $namepath . '.control_fragment.bw';
			}
		}
		
		if ($control_scale) {
			$self->{control_scale} = [split(',', $control_scale)];
			$self->crash("unequal scale factors and bam files!\n") if 
				scalar(@bams) != scalar(@{$self->{control_scale}});
		}
		# generate count bw  and dedup file names
		foreach my $bam (@bams) {
			my (undef, undef, $fname) = File::Spec->splitpath($bam);
			$fname =~ s/\.bam$//i; # strip extension
			my $base = File::Spec->catfile($opts{dir}, $fname);
			# count bigWig
			push @{ $self->{control_count_bw} }, "$base.count.bw";
			# dedup bam
			push @{ $self->{control_dedup_bams} }, "$base.dedup.bam";
		}
	}
	else {
		# must be just a chip without corresponding control
		$self->{control_bams} = [];
		$self->{lambda_bdg} = "$namepath.fragment.global_mean.bdg";
		$self->{lambda_bw} = "$namepath.fragment.global_mean.bw";
	}
	
	return $self;
}

sub crash {
	my ($self, $message) = @_;
	printf STDERR "ChIP Job object:\n%s\n", Dumper($self);
	confess $message;
}


sub generate_dedup_commands {
	my $self = shift;
	my $name2done = shift;
	croak("no bam_partial_dedup.pl script in path!\n") unless $opts{bamdedup} =~ /\w+/;
	my @commands;
	if (defined $self->{chip_bams}) {
		for (my $i = 0; $i < scalar @{$self->{chip_bams}}; $i++) {
			my $in = $self->{chip_bams}->[$i];
			my $out = $self->{chip_dedup_bams}->[$i];
			my $command = sprintf "%s --in %s --out %s --cpu %s ", 
				$opts{bamdedup},
				$in,
				$out,
				$opts{cpu};
			if (defined $opts{maxdepth} and $opts{maxdepth} == 1) {
				# no other deduplication options need to be set
				$command .= "--max 1 ";
			}
			else {
				# set random deduplication and/or maximum duplicates
				if ($opts{dupfrac} > 0) {
					$command .= sprintf("--seed 1 --frac %s ", $opts{dupfrac});
				}
				if (defined $opts{maxdepth} and $opts{maxdepth} > 1) {
					$command .= sprintf("--max %s ", $opts{maxdepth});
				}
			}
			if ($opts{optdist}) {
				$command .= sprintf("--optical --distance %s ", $opts{optdist});
			}
			if ($opts{paired} or $opts{deduppair}) {
				$command .= "--pe ";
			}
			if ($opts{blacklist}) {
				$command .= sprintf("--blacklist %s ", $opts{blacklist});
			}
			if ($opts{chrskip}) {
				$command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
			}
			my $log = $out;
			$log =~ s/\.bam$/.out.txt/i;
			$command .= " 2>&1 > $log";
			push @commands, [$command, $out, $log];
		}
	}
	if (defined $self->{control_bams}) {
		for (my $i = 0; $i < scalar @{$self->{control_bams}}; $i++) {
			my $in = $self->{control_bams}->[$i];
			if (exists $name2done->{$in}) {
				# this file has already been done, but we need to update the name
				$self->{control_dedup_bams}->[$i] = $name2done->{$in};
				next;
			}
			my $out = $self->{control_dedup_bams}->[$i];
			my $command = sprintf "%s --in %s --out %s --cpu %s ", 
				$opts{bamdedup},
				$in,
				$out,
				$opts{cpu};
			if (defined $opts{maxdepth} and $opts{maxdepth} == 1) {
				# no other deduplication options need to be set
				$command .= "--max 1 ";
			}
			else {
				# set random deduplication and/or maximum duplicates
				if ($opts{dupfrac} > 0) {
					$command .= sprintf("--seed 1 --frac %s ", $opts{dupfrac});
				}
				if (defined $opts{maxdepth} and $opts{maxdepth} > 1) {
					$command .= sprintf("--max %s ", $opts{maxdepth});
				}
			}
			if ($opts{optdist}) {
				$command .= sprintf("--optical --distance %s ", $opts{optdist});
			}
			if ($opts{paired} or $opts{deduppair}) {
				$command .= "--pe ";
			}
			if ($opts{blacklist}) {
				$command .= sprintf("--blacklist %s ", $opts{blacklist});
			}
			if ($opts{chrskip}) {
				$command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
			}
			my $log = $out;
			$log =~ s/\.bam$/.out.txt/i;
			$command .= " 2>&1 > $log";
			push @commands, [$command, $out, $log];
			$name2done->{$in} = $out; # remember that this has been done
		}
	}
	return @commands;
}

sub find_dedup_bams {
	my $self = shift;
	
	# ChIP bams
	if (defined $self->{chip_bams}) {
		printf " Checking ChIP bams for %s:\n", $self->{name};
		for (my $i = 0; $i < scalar @{$self->{chip_bams}}; $i++) {
			my $in = $self->{chip_bams}->[$i];
			my $out = $self->{chip_dedup_bams}->[$i];
			if (-e $out and -s _) {
				push @{$self->{chip_use_bams}}, $out;
				printf "  Found $out\n";
			}
			else {
				push @{$self->{chip_use_bams}}, $in;
				printf "  Using $in\n";
			}
		}	
	}
	
	# Control bams
	if (defined $self->{control_bams}) {
		printf " Checking control bams for %s:\n", $self->{name};
		for (my $i = 0; $i < scalar @{$self->{control_bams}}; $i++) {
			my $in = $self->{control_bams}->[$i];
			my $out = $self->{control_dedup_bams}->[$i];
			if (-e $out and -s _) {
				push @{$self->{control_use_bams}}, $out;
				printf "  Found $out\n";
			}
			else {
				push @{$self->{control_use_bams}}, $in;
				printf "  Using $in\n";
			}
		}
	}	
}

sub generate_bam2wig_frag_commands {
	my $self = shift;
	my $name2done = shift;
	croak("no bam2wig.pl script in path!\n") unless $opts{bam2wig} =~ /\w+/;
	my @commands;
	
	### ChIP bams
	if (scalar @{$self->{chip_use_bams}}) {
		# we have bam files to convert to bw
		
		# base fragment command
		my $frag_command = sprintf(
			"%s --out %s --qual %s --nosecondary --noduplicate --nosupplementary --bin %s --cpu %s --mean --bw --bwapp %s ", 
			$opts{bam2wig}, 
			$self->{chip_bw},
			$opts{mapq},
			$opts{chipbin},
			$opts{cpu},
			$opts{wig2bw}
		);
		
		# paired or single options
		if ($opts{paired}) {
			$frag_command .= sprintf("--pe --span --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
		}
		else {
			$frag_command .= sprintf("--extend --extval %s ", $opts{fragsize});
			$frag_command .= sprintf("--shiftval %s ", $opts{shiftsize}) 
				if $opts{shiftsize};
		}
		# additional filtering
		if ($opts{fraction}) {
			$frag_command .= "--fraction ";
		}
		if ($opts{blacklist}) {
			$frag_command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$frag_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		# scaling
		if ($self->{chip_scale}) {
			$frag_command .= sprintf("--scale %s ", join(',', @{$self->{chip_scale}} ) );
		}
		else {
			# standard scaling
			$frag_command .= "--rpm ";
		}
		# chromosome-specific scaling
		if ($self->{chrnorm}) {
			$frag_command .= sprintf("--chrnorm %s --chrapply %s ", $self->{chrnorm}, 
				$opts{chrapply});
		}
		# finish fragment command
		$frag_command .= sprintf("--in %s ", join(',', @{$self->{chip_use_bams}}));
		my $log = $self->{chip_bw};
		$log =~ s/bw$/out.txt/;
		$frag_command .= " 2>&1 > $log";
		push @commands, [$frag_command, $self->{chip_bw}, $log];
		
	}
	
	### Control bams
	my $control_bam_string = join(',', @{$self->{control_use_bams}});
	
	# first check if we've processed this control yet
	if (exists $name2done->{$control_bam_string}) {
		# we have, so change the lambda bedgraph name to the file that's been used
		# this will look funny, because it will look like we're using another chip's 
		# control file, but that's the way it will be
		$self->{lambda_bdg} = $name2done->{$control_bam_string};
	}
	
	# process control bams into a chromatin bias lambda-control track
	if (scalar @{$self->{control_use_bams}} and $opts{lambda} and 
		not exists $name2done->{$control_bam_string}
	) {
		
		# base fragment command
		my $frag_command = sprintf(
			"--qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --mean --bdg ", 
			$opts{mapq},
			$opts{cpu}
		);
		
		# general paired options, restrict size for all
		if ($opts{paired}) {
			$frag_command .= sprintf("--pe --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
		}
		
		# additional filters
		if ($opts{fraction}) {
			$frag_command .= "--fraction ";
		}
		if ($opts{blacklist}) {
			$frag_command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$frag_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		# chromosome-specific scaling
		if ($self->{chrnorm}) {
			$frag_command .= sprintf("--chrnorm %s --chrapply %s ", $self->{chrnorm}, 
				$opts{chrapply});
		}
		
		## now duplicate the base command for each size lambda control
		
		# d control, use extend or paired span just like ChIP
		my $command1 = $frag_command;
		if ($opts{paired}) {
			$command1 .= "--span ";
		}
		else {
			# treat d just like ChIP, including extend and shift, this may deviate from Macs2
			$command1 .= sprintf("--extend --extval %s --shiftval %0.0f --bin %s ", 
				$opts{fragsize}, $opts{shiftsize}, $opts{chipbin});
		}
		if ($self->{control_scale}) {
			$command1 .= sprintf("--scale %s ", join(',', @{$self->{control_scale}} ) );
		}
		else {
			$command1 .= "--rpm ";
		}
		# regenerate command 1 with program and input/output files
		$command1 = sprintf("%s --out %s %s --bin %s --in %s ", 
			$opts{bam2wig}, 
			$self->{d_control_bdg}, 
			$command1, 
			$opts{chipbin}, 
			$control_bam_string
		);
		my $log = $self->{d_control_bdg};
		$log =~ s/bdg$/out.txt/;
		$command1 .= " 2>&1 > $log";
		push @commands, [$command1, $self->{d_control_bdg}, $log];
		
		# small local lambda, extend both directions, scaled to compensate for length
		if ($opts{slocal}) {
			my $command2 = $frag_command;
			my $scale = sprintf("%.4f", $opts{fragsize}/$opts{slocal});
			if ($self->{control_scale}) {
				# user provided scale, multiply this with the lambda size scale
				my @scales = map {sprintf("%.4f", $_ * $scale)} @{$self->{control_scale}};
				$scale = join(',', @scales);
			}
			else {
				# standard scaling
				$command2 .= "--rpm "; 
			}
			# regenerate command 2 with program and input/output files
			$command2 = sprintf("%s --out %s %s --cspan --extval %s --scale %s --bin %s --in %s ", 
				$opts{bam2wig}, 
				$self->{s_control_bdg}, 
				$command2, 
				$opts{slocal}, 
				$scale, 
				$opts{slocalbin}, 
				$control_bam_string
			);
			$log = $self->{s_control_bdg};
			$log =~ s/bdg$/out.txt/;
			$command2 .= " 2>&1 > $log";
			push @commands, [$command2, $self->{s_control_bdg}, $log];
		}
		# large local lambda, extend both directions, scaled to compensate for length
		if ($opts{llocal}) {
			my $command3 = $frag_command;
			my $scale = sprintf("%.4f", $opts{fragsize}/$opts{llocal});
			if ($self->{control_scale}) {
				# user provided scale, multiply this with the lambda size scale
				my @scales = map {sprintf("%.4f", $_ * $scale)} @{$self->{control_scale}};
				$scale = join(',', @scales);
			}
			else {
				# standard scaling
				$command3 .= "--rpm "; 
			}
			# regenerate command 3 with program and input/output files
			$command3 = sprintf("%s --out %s %s --cspan --extval %s --scale %s --bin %s --in %s ", 
				$opts{bam2wig}, 
				$self->{l_control_bdg}, 
				$command3, 
				$opts{llocal}, 
				$scale, 
				$opts{llocalbin}, 
				$control_bam_string
			);
			$log = $self->{l_control_bdg};
			$log =~ s/bdg$/out.txt/;
			$command3 .= " 2>&1 > $log";
			push @commands, [$command3, $self->{l_control_bdg}, $log];
		}
		
		# record that we've done this bam
		# we store the lambda bedgraph name because it might be reused again for another job
		$name2done->{$control_bam_string} = $self->{lambda_bdg};
	}
	
	
	# skipping chromatin-bias lambda control track, use control track as is
	elsif (scalar @{$self->{control_use_bams}} and not $opts{lambda} and 
			not exists $name2done->{$control_bam_string}
	) {
		
		# fragment command
		my $frag_command = sprintf(
			"%s --out %s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --mean --bdg ", 
			$opts{bam2wig}, 
			$self->{lambda_bdg},
			$opts{mapq},
			$opts{cpu}
		);
		
		# single or paired options
		if ($opts{paired}) {
			$frag_command .= sprintf("--span --pe --minsize %s --maxsize %s --bin %s ", 
				$opts{minsize}, $opts{maxsize}, $opts{chipbin});
		}
		else {
			$frag_command .= sprintf("--extend --extval %s --shiftval %0.0f --bin %s ", 
				$opts{fragsize}, $opts{shiftsize}, $opts{chipbin});
		}
		
		# scaling for the fragment command only
		if ($self->{control_scale}) {
			$frag_command .= sprintf("--scale %s ", join(',', @{$self->{control_scale}} ) );
		}
		else {
			# standard scaling
			$frag_command .= "--rpm ";
		}
		
		# additional filters
		if ($opts{blacklist}) {
			$frag_command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$frag_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		# chromosome-specific scaling
		if ($self->{chrnorm}) {
			$frag_command .= sprintf("--chrnorm %s --chrapply %s ", $self->{chrnorm}, 
				$opts{chrapply});
		}
		# finish command
		$frag_command .= "--in $control_bam_string ";
		my $log = $self->{lambda_bdg};
		$log =~ s/bdg$/out.txt/;
		$frag_command .= " 2>&1 > $log";
		push @commands, [$frag_command, $self->{lambda_bdg}, $log];
		
		# record that we've done this bam
		# we store the lambda bedgraph name because it might be reused again for another job
		$name2done->{$control_bam_string} = $self->{lambda_bdg};
	}
	
	# finished
	return @commands;
}

sub generate_bam2wig_count_commands {
	my $self = shift;
	my $name2done = shift;
	croak("no bam2wig.pl script in path!\n") unless $opts{bam2wig} =~ /\w+/;
	my @commands;
	
	### ChIP bams
	if (scalar @{$self->{chip_use_bams}}) {
		# we have bam files to count
		
		# base count command
		my $count_command = sprintf(
			"--qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --bw --bwapp %s ", 
			$opts{mapq},
			$opts{cpu},
			$opts{wig2bw}
		);
		# paired or single options
		if ($opts{paired}) {
			$count_command .= sprintf("--pe --mid --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
		}
		else {
			$count_command .= sprintf("--start --shiftval %0.0f ", 
				($opts{fragsize} / 2) + $opts{shiftsize} );
		}
		# additional filtering
		if ($opts{fraction}) {
			$count_command .= "--fraction ";
		}
		if ($opts{blacklist}) {
			$count_command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$count_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		
		# finish count commands
		for my $i (0 .. $#{$self->{chip_use_bams}} ) {
			my $command = $count_command;
			# add scaling as necessary
			unless ($opts{rawcounts}) {
				if ($self->{chip_scale}) {
					$command .= sprintf("--scale %s ", $self->{chip_scale}->[$i]);
				}
				else {
					# we scale the count back up to target depth
					$command .= sprintf("--rpm --scale %d ", $opts{targetdep})
				}
				if ($self->{chrnorm}) {
					# chromosome-specific scaling
					$count_command .= sprintf("--chrnorm %s --chrapply %s ", 
						$self->{chrnorm}, $opts{chrapply});
				}
				$count_command .= "--format 3 "; # always limit digits
			}
			# regenerate command with file names
			$command = sprintf("%s --out %s %s --in %s ", 	
				$opts{bam2wig}, 
				$self->{chip_count_bw}->[$i], 
				$command,
				$self->{chip_use_bams}->[$i]
			);
			# output log
			my $log = $self->{chip_count_bw}->[$i];
			$log =~ s/bw$/out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [$command, $self->{chip_count_bw}->[$i], $log];
		}
	}
	
	### Control bams
	# we only process this if we haven't done them yet
	my $control_bam_string = join(',', @{$self->{control_use_bams}});
	if (scalar @{$self->{control_use_bams}} and 
		not exists $name2done->{$control_bam_string}
	) {
		
		# base count command
		my $count_command = sprintf(
			"--qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --bw --bwapp %s ", 
			$opts{mapq},
			$opts{cpu},
			$opts{wig2bw}
		);
		
		# general paired options, restrict size for all
		if ($opts{paired}) {
			$count_command .= sprintf("--pe --mid --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
		}
		else {
			$count_command .= sprintf("--start --shiftval %0.0f ", 
				($opts{fragsize} / 2) + $opts{shiftsize});
		}
		
		# additional filters
		if ($opts{fraction}) {
			$count_command .= "--fraction ";
		}
		if ($opts{blacklist}) {
			$count_command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$count_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		
		# add count commands
		for my $i (0 .. $#{$self->{control_use_bams}} ) {
			my $command = $count_command;
			# add scaling as necessary
			unless ($opts{rawcounts}) {
				if ($self->{control_scale}) {
					$command .= sprintf("--scale %s ", $self->{control_scale}->[$i]);
				}
				else {
					# we scale the count back up to target depth
					$command .= sprintf("--rpm --scale %d ", $opts{targetdep})
				}
				if ($self->{chrnorm}) {
					# chromosome-specific scaling
					$count_command .= sprintf("--chrnorm %s --chrapply %s ", 
						$self->{chrnorm}, $opts{chrapply});
				}
				$command .= "--format 3 ";
			}
			# regenerate command with file names
			$command = sprintf("%s --out %s %s --in %s ", 	
				$opts{bam2wig}, 
				$self->{control_count_bw}->[$i], 
				$command,
				$self->{control_use_bams}->[$i]
			);
			# output log
			my $log = $self->{control_count_bw}->[$i];
			$log =~ s/bw$/out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [$command, $self->{control_count_bw}->[$i], $log];
		}
		
		# record that we've done this bam
		$name2done->{$control_bam_string} = 1;
	}
	
	# finished
	return @commands;
}


sub generate_lambda_control_commands {
	my $self = shift;
	my $name2done = shift;
	return if exists $name2done->{ $self->{lambda_bdg} }; # already done
	
	# check whether no reference control was provided at all
	if (not scalar @{$self->{control_use_bams}} and 
		$self->{lambda_bdg} =~ /global_mean\.bdg$/ and
		not -e $self->{lambda_bdg}
	) {
		# no control bam files at all, need to use expected mean
		# generate this from the existing ChIP fragment bigWig file
		# no output file needed to be specified, it will automatically append .global_mean.bdg
		my $log = $self->{lambda_bdg};
		$log =~ s/bdg/out.txt/;
		my $command = sprintf("%s %s 2>&1 > $log", $opts{meanbdg}, $self->{chip_bw});
		$name2done->{ $self->{lambda_bdg} } = 1;
		return [$command, $self->{lambda_bdg}, $log];
	}
	
	# Proceed with generating lambda bedGraph
	croak("no macs2 application in path!\n") unless $opts{macs} =~ /\w+/;
	croak("no bedtools application in path!\n") unless $opts{bedtools} =~ /\w+/;
	my $dfile = $self->{d_control_bdg};
	my $sfile = $self->{s_control_bdg};
	my $lfile = $self->{l_control_bdg};
	return unless ($dfile); # controls with lambda will always have  d_control_bdg
		# ChIP jobs and non-lambda controls will not
	unless ($opts{dryrun}) {
		$self->crash("no d control bedGraph file '$dfile'!\n") if (not -e $dfile);
		$self->crash("no small control bedGraph file '$sfile'!\n") if ($sfile and not -e $sfile);
		$self->crash("no large control bedGraph file '$lfile'!\n") if ($lfile and not -e $lfile);
	}
	
	# generate background lambda bedgraph
	# go ahead and do this immediately because it's so simple
	my $background_bdg;
	unless ($opts{dryrun}) {
		my $background = sprintf("%.4f", ( 1_000_000 * $opts{fragsize} ) / $opts{genome} );
		printf " Calculating background for %s as $background\n", $self->{name};
		$background_bdg = $self->{lambda_bdg};
		$background_bdg =~ s/lambda_control/background/;
		my $infh = IO::File->new($chromofile, 'r') or  # use the chromosome file as source
			die "unable to open chromosome file '$chromofile'!\n";
		my $outfh = IO::File->new($background_bdg, "w");
		while (my $line = $infh->getline) {
			chomp $line;
			my ($chr, $end) = split /\s/, $line;
			$outfh->printf("%s\t0\t%s\t%s\n", $chr, $end, $background);
		}
		$infh->close;
		$outfh->close;
	}
	
	my $log = $self->{lambda_bdg};
	$log =~ s/bdg$/out.txt/;
	my $command;
	
	# generate commands using bedtools and data2wig, which is faster than running 
	# multiple instances of macs2 bdgcmp and bdgopt
	if ($sfile and $lfile) {
		# first step
		$command = sprintf("%s unionbedg -header -names dlocal slocal llocal background -i %s %s %s %s > %s 2> $log ", 
			$opts{bedtools}, $dfile, $sfile, $lfile, $background_bdg, 
			$self->{sld_control_file});
		
		# second step
		$command .= sprintf("&& %s --in %s --zero --fast --bdg --notrack --score 3-6 --method max --out %s ",
			$opts{data2wig}, $self->{sld_control_file}, $self->{lambda_bdg});
		$command .= " 2>&1 >> $log ";
		
		# clean up
		$command .= sprintf("&& rm %s %s %s %s ", $sfile, $lfile, 
			$background_bdg, $self->{sld_control_file});
	}
	elsif ($sfile and not $lfile) {
		# first step
		$command = sprintf("%s unionbedg -header -names dlocal slocal background -i %s %s %s > %s 2> $log ", 
			$opts{bedtools}, $dfile, $sfile, $background_bdg, 
			$self->{sld_control_file});
		
		# second step
		$command .= sprintf("&& %s --in %s --zero --fast --bdg --notrack --score 3-5 --method max --out %s ",
			$opts{data2wig}, $self->{sld_control_file}, $self->{lambda_bdg});
		$command .= " 2>&1 >> $log ";
		
		# clean up
		$command .= sprintf("&& rm %s %s %s ", $sfile, 
			$background_bdg, $self->{sld_control_file});
	}
	elsif (not $sfile and $lfile) {
		# first step
		$command = sprintf("%s unionbedg -header -names dlocal llocal background -i %s %s %s > %s 2> $log ", 
			$opts{bedtools}, $dfile, $lfile, $background_bdg, 
			$self->{sld_control_file});
		
		# second step
		$command .= sprintf("&& %s --in %s --zero --fast --bdg --notrack --score 3-5 --method max --out %s ",
			$opts{data2wig}, $self->{sld_control_file}, $self->{lambda_bdg});
		$command .= " 2>&1 >> $log ";
		
		# clean up
		$command .= sprintf("&& rm %s %s %s ", $lfile, 
			$background_bdg, $self->{sld_control_file});
	}
	else {
		$self->crash("programming error! how did we get here with no sfile and no lfile????");
	}
	
	$name2done->{ $self->{lambda_bdg} } = 1;
	return [$command, $self->{lambda_bdg}, $log];
}

sub convert_bw_to_bdg {
	my $self = shift;
	my $name2done = shift;
	croak("no bigWigToBedGraph application in path!\n") unless $opts{bw2bdg} =~ /\w+/;
	my @commands;
	if (
		($self->{chip_bw} and -e $self->{chip_bw}) or 
		$opts{dryrun}
	) {
		# generate command if bigWig file exists or we're in dry run mode
		my $log = $self->{chip_bdg};
		$log =~ s/bdg$/out.txt/;
		my $command = sprintf("%s %s %s 2>> $log", $opts{bw2bdg}, $self->{chip_bw}, 
			$self->{chip_bdg});
		push @commands, [$command, $self->{chip_bdg}, $log]
	}
	if (
		(
			$self->{lambda_bw} and -e $self->{lambda_bw} and 
			not exists $name2done->{$self->{lambda_bdg}} 
		) or
		$opts{dryrun}
	) {
		# generate command if bigWig exists and hasn't been done yet
		# or if we're in dry run mode
		my $log = $self->{lambda_bdg};
		$log =~ s/bdg$/out.txt/;
		my $command = sprintf("%s %s %s 2>> $log", $opts{bw2bdg}, $self->{lambda_bw}, 
			$self->{lambda_bdg});
		push @commands, [$command, $self->{lambda_bdg}, $log];
		$name2done->{$self->{lambda_bdg}} = 1; # finished
	}
	return @commands;
}

sub generate_enrichment_commands {
	my $self = shift;
	croak("no macs2 application in path!\n") unless $opts{macs} =~ /\w+/;
	my $chip = $self->{chip_bdg} || undef;
	my $lambda = $self->{lambda_bdg} || undef;
	return unless ($chip and $lambda);
	unless ($opts{dryrun}) {
		$self->crash("no ChIP fragment file $chip!\n") if not -e $chip;
		$self->crash("no control lambda fragment file $lambda!\n") if not -e $lambda;
	}
	
	my $command = sprintf("%s bdgcmp -t %s -c %s -m qpois FE ", $opts{macs}, $chip, 
		$lambda);
	if (defined $opts{targetdep}) {
		$command .= sprintf("-S %d ", $opts{targetdep});
	}
	if (not $opts{lambda}) {
		$command .= "-p 1 "; # add a pseudo count of 1 when doing explicit comparisons
	}
	$command .= sprintf("-o %s %s ", $self->{qvalue_bdg}, $self->{fe_bdg});
	my $log = $self->{qvalue_bdg};
	$log =~ s/qvalue\.bdg$/bdgcmp.out.txt/;
	$command .= " 2> $log ";
	return [$command, $self->{qvalue_bdg}, $log];
}

sub generate_peakcall_commands {
	my $self = shift;
	croak("no macs2 application in path!\n") unless $opts{macs} =~ /\w+/;
	my $qtrack = $self->{qvalue_bdg} || undef;
	return unless $qtrack;
	unless ($opts{dryrun}) {
		$self->crash("no qvalue bedGraph file $qtrack!\n") if not -e $qtrack;
	}
	my $command = sprintf("%s bdgpeakcall -i %s -c %s -l %s -g %s --no-trackline -o %s ",
		$opts{macs}, $qtrack, $opts{cutoff}, $opts{peaksize}, $opts{peakgap}, 
		$self->{peak}
	);
	my $log = $self->{peak};
	$log =~ s/narrowPeak$/peakcall.out.txt/;
	$command .= " 2> $log";
	if ($opts{broad}) {
		# sneak this in as an extra command
		# bdgbroadcall doesn't support writing without a trackline
		my $command2 = sprintf("%s bdgbroadcall -i %s -c %s -C %s -l %s -g %s -G %s -o %s ",
			$opts{macs}, $qtrack, $opts{cutoff}, $opts{broadcut}, $opts{peaksize}, 
			$opts{peakgap}, $opts{broadgap}, $self->{gappeak}
		);
		my $log2 = $self->{gappeak};
		$log2 =~ s/gappedPeak$/broadcall.out.txt/;
		$command2 .= " 2> $log2";
		return ( [$command, $self->{peak}, $log], [$command2, $self->{gappeak}, $log2] );
	}
	return [$command, $self->{peak}, $log];
}

sub generate_cleanpeak_commands {
	my $self = shift;
	return unless defined $self->{peak}; # skip control jobs
	croak("no peak2bed.pl application in path!\n") unless $opts{peak2bed} =~ /\w+/;
	
	# this is now done by an external script
	# we do not specify output file names, but it uses the same logic as here
	# there are no options to this script
	my $command1 = sprintf("%s ", $opts{peak2bed});
	my $command2 = 'rm ';
	my $command_check = length($command1);
	
	# narrowPeak
	if (
		(-e $self->{peak} and -s _ > 0) or 
		$opts{dryrun}
	) {
		# generate the command only if the peak file exists and nonzero in size
		# or we're in dry run mode
		$command1 .= sprintf("%s ", $self->{peak});
		$command2 .= sprintf("%s ", $self->{peak});
	}
	# gappedPeak
	if (
		$opts{broad} and 
		(-e $self->{gappeak} or $opts{dryrun})
	) {
		# generate the command only if the gapped peak file exists
		# or we're in dry run mode
		$command1 .= sprintf("%s ", $self->{gappeak});
		$command2 .= sprintf("%s ", $self->{gappeak});
	}
	
	# check that we have added to the command
	if (length($command1) > $command_check) {
		# good, we have output files and are doing something
		my $log = $self->{clean_peak};
		$log =~ s/bed/cleanpeak.out.txt/;
		$command1 .= sprintf(" 2>&1 > %s ", $log);
		$command1 .= sprintf("&& %s", $command2);
		return [$command1, $self->{clean_peak}, $log];
	}
	else {
		# empty file(s), let's fake the clean one
		if ($opts{broad}) {
			$command1 = sprintf "touch %s %s %s && rm -f %s %s", $self->{clean_peak}, 
				$self->{peak_summit}, $self->{clean_gappeak}, $self->{peak}, 
				$self->{gappeak};
		}
		else {
			$command1 = sprintf "touch %s %s && rm %s", $self->{clean_peak}, 
				$self->{peak_summit}, $self->{peak};
		}
		return [$command1, $self->{clean_peak}, ''];
	}
}

sub generate_bdg2bw_commands {
	my $self = shift;
	my $chromofile = shift;
	my $name2done = shift;
	croak("no wigToBigWig application in path!\n") unless $opts{wig2bw} =~ /\w+/;
	croak("no manipulate_wig.pl script in path!\n") unless $opts{manwig} =~ /\w+/;
	
	my @commands;
	if ($self->{chip_bdg} and $self->{chip_bw}) {
		# we should have both files here
		if (-e $self->{chip_bdg} and -e $self->{chip_bw}) {
			# bigWig already exists so just delete the bdg
			push @commands, [sprintf("rm %s",$self->{chip_bdg}), '', ''];
		}
	}
	if ($self->{lambda_bdg} and $self->{lambda_bw} and 
		not exists $name2done->{$self->{lambda_bdg}}
	) {
		if (-e $self->{lambda_bw} ) {
			# we must have started with a lambda bigWig so remove the bedGraph
			push @commands, [sprintf("rm %s", $self->{lambda_bdg}), '', ''];
			$name2done->{$self->{lambda_bdg}} = 1;
		}
		else {
			# lambda control bigWig
			my $log = $self->{lambda_bw};
			$log =~ s/bw$/out.txt/;
			my $command = sprintf("%s %s %s %s 2>&1 > $log && rm %s", 
				$opts{wig2bw},
				$self->{lambda_bdg},
				$chromofile,
				$self->{lambda_bw},
				$self->{lambda_bdg},
			);
			push @commands, [$command, $self->{lambda_bw}, $log];
			
			# d control fragment bigWig
			if ($self->{d_control_bdg} and $self->{d_control_bw}) {
				$log = $self->{d_control_bw};
				$log =~ s/bw$/out.txt/;
				$command  = sprintf("%s %s %s %s 2>&1 > $log && rm %s", 
					$opts{wig2bw},
					$self->{d_control_bdg},
					$chromofile,
					$self->{d_control_bw},
					$self->{d_control_bdg},
				);
				push @commands, [$command, $self->{d_control_bw}, $log];
			}
		}
		$name2done->{$self->{lambda_bdg}} = 1; # remember it's done
	}
	if ($self->{qvalue_bdg} and $self->{qvalue_bw}) {
		my $log = $self->{qvalue_bw};
		$log =~ s/bw$/out.txt/;
		my $command = sprintf("%s %s %s %s 2>&1 > $log ", 
			$opts{wig2bw},
			$self->{qvalue_bdg},
			$chromofile,
			$self->{qvalue_bw}
		);
		unless ($opts{savebdg}) {
			$command .= sprintf("&& rm %s", $self->{qvalue_bdg});
		}
		push @commands, [$command, $self->{qvalue_bw}, $log];
	}
	if ($self->{fe_bdg}) {
		# convert this to log2 Fold Enrichment because I like this better
		my $log = $self->{logfe_bw};
		$log =~ s/bw$/out.txt/;
		my $command = sprintf("%s --in %s --log 2 --place 4 --w2bw %s --chromo %s --out %s 2>&1 > $log ",
			$opts{manwig},
			$self->{fe_bdg},
			$opts{wig2bw},
			$chromofile,
			$self->{logfe_bw},
		);
		$command .= sprintf(" && rm %s", $self->{fe_bdg});
		push @commands, [$command, $self->{logfe_bw}, $log];
	}
	return @commands;
}





