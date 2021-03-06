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
use File::Spec;
use File::Path qw(make_path);
use List::Util qw(uniqstr);
use Getopt::Long;
use Bio::MultiRepChIPSeq::Runner;

# Initialize Runner and options
my $Runner = Bio::MultiRepChIPSeq::Runner->new();
my $opts = $Runner->options;
my $VERSION = $Runner->version;

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
fragment (or d in Macs2 parlance), small lambda (default $opts->{slocal} bp), and 
large lambda (default $opts->{llocal} bp) fragment coverage. If desired, either small or 
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
  --dir       directory         Directory for writing all files ($opts->{dir})
  --out       file basename     Base filename for merged output files ($opts->{out})
  
 Genome size
  --genome    integer           Specify effective mappable genome size 
                                  (default empirically determined)
  
 Bam options
  --mapq      integer           Minimum mapping quality, ($opts->{mapq})
  --pe                          Bam files are paired-end, default treat as single-end
  --min       integer           Minimum paired-end size allowed ($opts->{minsize} bp)
  --max       integer           Maximum paired-end size allowed ($opts->{maxsize} bp)
  --fraction                    Record multiple-hit alignments as fraction of hits
 
 Bam filtering options
  --chrskip   "text"            Chromosome skip regex ($opts->{chrskip})
  --blacklist file              Bed file of repeats or hotspots to avoid
                                  Determined empirically from control (Input) samples
  
 Duplication filtering
  --nodedup                     Skip deduplication and take everything as is
  --dupfrac   float             Target duplication rate for subsampling ($opts->{dupfrac})
  --maxdepth  integer           Maximum position alignment depth ($opts->{maxdepth})
                                  set to 1 to remove all duplicates
  --optdist   integer           Maximum distance for optical duplicates ($opts->{optdist})
                                  use 100 for HiSeq, 2500 for NovaSeq
  --deduppair                   Run deduplication as paired-end, but coverage as single-end
                                  e.g. for ATAC-Seq cut site analysis

 Fragment coverage
  --size      integer           Predicted fragment size. REQUIRED for single-end
  --shift     integer           Shift the fragment, e.g. ATACSeq ($opts->{shiftsize} bp)
  --slocal    integer           Small local lambda size ($opts->{slocal} bp)
  --llocal    integer           Large local lambda size ($opts->{llocal} bp)
  --cbin      integer           ChIP fragment bin size ($opts->{chipbin} bp)
  --slbin     integer           Small local lambda bin size ($opts->{slocalbin} bp)
  --llbin     integer           Large local lambda bin size ($opts->{llocalbin} bp)

 Chromosome-specific normalization
  --chrnorm   float             Specific chromosome normalization factor
  --chrapply  "text"            Apply factor to specified chromosomes via regex
 
 Peak calling
  --cutoff    number            Threshold q-value for calling peaks ($opts->{cutoff}) 
                                  Higher numbers are more significant, -1*log10(q)
  --peaksize  integer           Minimum peak size to call (2 x fragment size)
                                  Required for paired-end alignments.
  --peakgap   integer           Maximum gap between peaks before merging (1 x size)
  --broad                       Also perform broad (gapped) peak calling
  --broadcut  number            Q-value cutoff for linking broad regions ($opts->{broadcut})
  --broadgap  integer           Maximum link size between peaks in broad calls (4 x size bp)
  --nolambda                    Skip lambda control, compare ChIP directly with control
  --independent                 Call peaks independently for each replicate and merge
  
 Peak scoring
  --binsize   integer           Size of bins in 25 flanking peak bins for profile ($opts->{binsize})
  --genomewin integer           Collect counts across genome in given window size
  --discard   number            Discard genome windows with replicate sum below number ($opts->{discard})
  --rawcounts                   Use unscaled raw counts for re-scoring peaks
  --repmean                     Combine replicate counts as mean for each sample set
  --noplot                      Do not plot figures of results
  
 Job control
  --cpu       integer           Number of CPUs to use per job ($opts->{cpu})
  --job       integer           Number of simultaneous jobs ($opts->{job})
  --dryrun                      Just print the commands without execution
  --noorganize                  Do not organize files into subfolders when finished
  --savebam                     Save de-duplicated bam files
  --savebdg                     Save text bedGraph files

 Application  Paths
  --bam2wig   path             ($opts->{bam2wig})
  --bamdedup  path             ($opts->{bamdedup})
  --macs      path             ($opts->{macs})
  --manwig    path             ($opts->{manwig})
  --wig2bw    path             ($opts->{wig2bw})
  --bw2bdg    path             ($opts->{bw2bdg})
  --printchr  path             ($opts->{printchr})
  --data2wig  path             ($opts->{data2wig})
  --getdata   path             ($opts->{getdata})
  --getrel    path             ($opts->{getrel})
  --geteff    path             ($opts->{geteff})
  --meanbdg   path             ($opts->{meanbdg})
  --bedtools  path             ($opts->{bedtools})
  --intersect path             ($opts->{intersect})
  --peak2bed  path             ($opts->{peak2bed})
  --combrep   path             ($opts->{combrep})
  --plotpeak  path             ($opts->{plotpeak})
  --rscript   path             ($opts->{rscript})
  --reportmap path             ($opts->{reportmap})
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
	$opts,
	'dir=s',
	'out=s',
	'chip=s@',
	'control=s@',
	'name=s@',
	'chscale=s@',
	'coscale=s@',
	'species=s',
	'genome=i',
	'mapq=i',
	'paired|pe!',
	'fraction!',
	'minsize|min=i',
	'maxsize|max=i',
	'chrskip=s',
	'blacklist=s',
	'dedup!',
	'dupfrac=f',
	'maxdup=i', # old option
	'maxdepth=i',
	'optdist=i',
	'deduppair!',
	'fragsize|size=i',
	'shiftsize|shift=i',
	'slocal=i',
	'llocal=i',
	'chipbin|cbin=i',
	'slocalbin|slbin=i',
	'llocalbin|llbin=i',
	'chrnorm=f@',
	'chrapply=s',
	'cutoff=f',
	'targetdep|tdep=f',
	'peaksize=i',
	'peakgap=i',
	'broad!',
	'broadcut=f',
	'broadgap=i',
	'lambda!',
	'independent!',
	'binsize=i',
	'genomewin=i',
	'discard=f',
	'repmean!',
	'plot!',
	'cpu=i',
	'job=i',
	'dryrun!',
	'organize!',
	'savebam!',
	'savebdg!',
	'help!',
	'bam2wig=s',
	'bamdedup=s',
	'macs=s',
	'manwig=s',
	'wig2bw=s',
	'bw2bdg=s',
	'printchr=s',
	'data2wig=s',
	'getdata=s',
	'getrel=s',
	'geteff=s',
	'meanbdg=s',
	'bedtools=s',
	'intersect=s',
	'peak2bed=s',
	'combrep=s',
	'plotpeak=s',
	'rscript=s',
	'reportmap=s',
) or die "unrecognized option(s)!\n";

if ($opts->{help}) {
	print $documentation;
	exit 1;
}



### Begin main pipeline
print "======== ChIPSeq multi-replicate pipeline ==========\n";
print "\nversion $VERSION\n";
my $start = time;
check_inputs();
print_start();
check_input_files();
add_jobs_to_runner();

# Check progress
$Runner->check_progress_file();

# Start pipeline
$Runner->run_generate_chr_file();
$Runner->run_input_peak_detection();
$Runner->run_dedup();
$Runner->run_bam_check();
$Runner->run_mappable_space_report();
$Runner->run_bam_fragment_conversion();
$Runner->run_bam_count_conversion();
$Runner->run_lambda_control();
$Runner->run_bw_conversion();
$Runner->run_bdgcmp();
$Runner->run_call_peaks();
$Runner->run_clean_peaks();
$Runner->run_peak_merge();
$Runner->run_bdg_conversion();
$Runner->write_samples_file();
$Runner->run_rescore();
$Runner->run_efficiency();
$Runner->run_plot_peaks();
$Runner->run_cleanup();
$Runner->run_organize();

# final statement
printf "\n\nFinished in %.1f minutes\n", (time -$start) / 60;
print "======== Finished ChIPSeq multi-replicate pipeline ==========\n";




############### Subroutines ########################################################

sub check_inputs {
	if (@ARGV) {
		die sprintf("There are unrecognized leftover items on the command line!\n Did you leave spaces in your --chip or --control file lists?\nItems:\n %s\n", join("\n ", @ARGV));
	}
	unless ($Runner->chip) {
		die "No ChIP file(s) defined!\n";
	}
	if ($Runner->name) {
		my %check = map {$_ => 1} ($Runner->name);
		if (scalar(keys %check) != scalar($Runner->name)) {
			die "Duplicate sample names are present!\n";
		}
	}
	else {
		die "No name(s) defined!\n";
	}
	unless (scalar($Runner->chip) == scalar($Runner->name)) {
		die "Unequal number of ChIP samples and names!\n";
	}
	
	if (scalar($Runner->control) > 1 and scalar($Runner->control) != scalar($Runner->chip)) {
		die "Unequal number of control and ChIP samples!\n";
	}
	elsif (scalar($Runner->control) == 0) {
		# no controls, turn off lambda
		$Runner->slocal(0);
		$Runner->llocal(0);
		$Runner->lambda(0);
	}
	if ($Runner->slocal == 0 and $Runner->llocal == 0) {
		# user somehow set both to zero, but this throws errors to macs2
		$Runner->lambda(0);
	}
	if (scalar($Runner->chscale) or scalar($Runner->coscale)) {
		# no longer recommended
			print <<MESSAGE;

WARNING: Manually setting ChIP and/or Control scaling factors!!!!
Please be aware that manually scaling coverage depth is an advanced option
and should only be done when you're aware of the ramifications. Inappropriate
scaling may artificially alter the expected statistics required for accurate
peak calling and may reduce confidence of identified peaks.

MESSAGE
	}
	if (scalar($Runner->chscale) and scalar($Runner->chscale) != scalar($Runner->chip)) {
		die "Unequal number of ChIP samples and ChIP scale factors!\n";
	}
	if (scalar($Runner->coscale) and scalar($Runner->coscale) != scalar($Runner->control)) {
		die "Unequal number of control samples and control scale factors!\n";
	}
	if (scalar($Runner->chrnorm) and not $Runner->chrapply) {
		die "Must specify chromosome apply regex (--chrapply) for chromosome normalization factors!\n";
	}
	if (scalar($Runner->chrnorm) and scalar($Runner->chrnorm) != scalar($Runner->chip)) {
		# apply to all the ChIPs
		if (scalar($Runner->chrnorm) == 1) {
			print "Using the same chromosome normalization factor for each ChIP sample\n";
			my $n = ($Runner->chrnorm)[0];
			my $m = scalar($Runner->name);
			for (2 .. $m) {
				# add it to the remainder
				$Runner->chrnorm($n);
			}
		}
		else {
			die "Unequal number of chromosome normalization factors and ChIP samples!\n";
		}
	}
	if (scalar($Runner->chrnorm) > 1 and scalar($Runner->control) == 1) {
		print "Using first chromosome normalization factor for universal control!\n";
	}
	if ($Runner->species) {
		print <<MESSAGE;

WARNING: Specifiying species is now deprecated. The genome mappable size is 
now determined empirically from all provided Bam files using the script 
report_mappable_space.pl, or an explicit genome mappable size may be provided 
with the --genome option. 

MESSAGE
	}
	
	# directory
	unless ($Runner->dryrun) {
		unless (-e $Runner->dir and -d _) {
			make_path($Runner->dir) or 
				die sprintf("Unable to make directory %s! $!\n", $Runner->dir);
		}
		unless (-w $Runner->dir) {
			die sprintf("Target directory %s is not writable!\n", $Runner->dir);
		}
	}
	# target depth
	if (defined $Runner->targetdep) {
		my $files = join(",", ($Runner->chip), ($Runner->control));
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
	if (not $Runner->peaksize) {
		# no minimum peak size defined? might be ok
		if ($Runner->fragsize) {
			# set default to twice fragment size
			$Runner->peaksize(2 * $Runner->fragsize);
		}
		elsif ($Runner->paired) {
			# paired fragments - make something up
			if (not $Runner->peaksize) {
				print "\nWARNING! Setting minimum peak size to 500 bp, but this should be manually\nset based on mean alignment insert size and nature of experiment.\n\n";
				$Runner->peaksize(500);
				$Runner->fragsize(250); # for expected background normalization
			}
		}
		else {
			die "Must set an estimated mean fragment size for single-end alignments!\n  Run 'macs2 predictd' or 'bam2wig.pl --shift --model'\n";
		}
	}
	unless (defined $Runner->peakgap) {
		my $g = $Runner->fragsize ? $Runner->fragsize : int($Runner->peaksize/2);
		$Runner->peakgap($g);
	}
	unless (defined $Runner->broadgap) {
		my $g = $Runner->fragsize ? 4 * $Runner->fragsize : 2 * $Runner->peaksize;
		$Runner->broadgap($g);
	}
	
	# exclusion list
	if (
		$Runner->blacklist and 
		$Runner->blacklist ne 'input' and 
		not -e $Runner->blacklist
	) {
		printf("\nWARNING! Unable to find specified black list file '%s'!\n", 
			$Runner->blacklist);
		if (scalar($Runner->control)) {
			print "Defaulting to using input-derived exclusion list\n";
			$Runner->blacklist('input');
		}
		else {
			$Runner->blacklist(q());
		}
	}
	if (not defined $Runner->blacklist and scalar($Runner->control)) {
		$Runner->blacklist('input');
	}
	
	# max depth-duplication confusion
	if ($Runner->maxdup) {
		# because this was inappropriately named before
		print " \nWARNING: The --maxdup option is now --maxdepth\n";
		$Runner->maxdepth($Runner->maxdup);
	}
}

sub print_start {
	if ($Runner->dryrun) {
		print <<DRYRUN;

======= Dry Run =======
Some values are empirically determined during execution and are made up here.
Some files may not be generated during actual execution, but commands should
mostly be complete. No files will be generated.
=======================

DRYRUN
	}
	$Runner->print_config;
	print "\n\n";
}

sub check_input_files {
	print "\n\n======= Checking input files\n";
	my @errors;
	my @names = $Runner->name;
	my @chips = $Runner->chip;
	my @controls = $Runner->control;
	foreach my $i (0 .. $#names) {
		my @list = split ',', $chips[$i];
		if (scalar(@list) != uniqstr(@list)) {
			printf " Fixing duplicate entries found in %s ChIP entries!\n", $names[$i];
			# fix it
			@list = uniqstr(@list);
			# hack to replace
			$chips[$i] = join(',', @list);
			$opts->{chip} = \@chips;
		}
		foreach my $f (@list) {
			unless (-e $f) {
				printf(" can't find %s ChIP file %s!\n", $names[$i], $f);
				push @errors, $f;
			}
		}
		
		if (@controls) {
			my @list2 = split ',', $controls[$i];
			if (scalar(@list2) != uniqstr(@list2)) {
				printf " Fixing duplicate entries found in %s control entries!\n", $names[$i];
				# fix it
				@list2 = uniqstr(@list2);
				# hack to replace
				$controls[$i] = join(',', @list2);
				$opts->{control} = \@controls;
			}
			foreach my $f (@list2) {
				unless (-e $f) {
					printf(" can't find %s control file %s!\n", $names[$i], $f);
					push @errors, $f;
				}
			}
		}
	}
	if (scalar @errors) {
		die sprintf("ERROR! Missing %d input files!\nMissing: %s\n", 
			scalar @errors, join(', ', @errors) );
	}
	else {
		print " All input files found\n";
	}
}

sub add_jobs_to_runner {
	# we pass this off to ChIPjob package with 6 options: 
	# name, chip files, control files, chip scale factor, control scale factor, 
	# chromosome normalization factor
	
	# first check for a unversal control
	if (scalar($Runner->control) == 1 and scalar($Runner->chip) > 1) {
		print "Using only one control for multiple ChIP experiments\n";
		my $universal_control = ($Runner->control)[0];
		
		# generate universal name
		my $universal_name;
		if ($universal_control =~ /\.bam$/i) {
			# one bam file
			$universal_name = $Runner->out . '_control';
			if ($Runner->lambda) {
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
		my $universal_scale = ($Runner->coscale)[0] || undef;
		
		# add back universal control with special prefix for each ChIP job
		$opts->{control} = []; # unorthodox hack
		$opts->{coscale} = []; # unorthodox hack
		foreach ($Runner->name) {
			$Runner->control('Custom-Universal-' . $universal_name);
		}
		
		# generate job
		$Runner->add_job(
			$universal_name, 
			'', 
			$universal_control, 
			undef, 
			$universal_scale, 
			($Runner->chrnorm)[0] || undef
		);
	}
	
	# walk through each given job
	my @names = $Runner->name;
	for my $i (0 .. $#names) {
		# check name
		my $name = $names[$i];
		if ($name =~ /\./) {
			# BioToolBox simplify names by dropping everything after a period
			# so a crude hack is to change the name
			$name =~ s/\./_/g; 
			printf "\nWARNING! Changing name from '%s' to '%s' for compatability\n", 
				$names[$i], $name;
			# ugly hack to replace the name
			$opts->{name}->[$i] = $name;
		}
		# add job
		$Runner->add_job(
			$name, 
			($Runner->chip)[$i], 
			($Runner->control)[$i] || undef, 
			($Runner->chscale)[$i] || undef, 
			($Runner->coscale)[$i] || undef, 
			($Runner->chrnorm)[$i] || undef
		);
	}
}



