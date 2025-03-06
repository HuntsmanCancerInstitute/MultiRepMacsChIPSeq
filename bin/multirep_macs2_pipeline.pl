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
use English qw(-no_match_vars);
use File::Spec;
use File::Path qw(make_path);
use List::Util qw(uniqstr);
use Getopt::Long;
use Bio::MultiRepChIPSeq::Runner;

# Initialize Runner and options
my $Runner = Bio::MultiRepChIPSeq::Runner->new();
my $opts   = $Runner->options;
our $VERSION = $Runner->version;

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
  --chip        file1,file2...  Repeat for each sample set
  --name        text            Repeat for each sample
  --control     file1,file2...  Repeat if matching multiple samples
 
 Output
  --dir         directory       Directory for writing all files ($opts->{dir})
  --out         file basename   Base filename for merged output files ($opts->{out})
  
 Genome size
  --genome      integer         Specify effective mappable genome size 
                                (default empirically determined)
  
 Bam options
  --pe                          Bam files are paired-end, default treat as single-end
 
 Alignment filtering options
  --mapq        integer         Minimum mapping quality, ($opts->{mapq})
  --chrskip     "text"          Chromosome skip regex ($opts->{chrskip})
  --exclude     file            Coordinate (bed) file of repeats or hotspots to avoid
                                  Default determined empirically from control samples.
                                  Specify 'none' for no filtering.
  --min         integer         Minimum paired-end size allowed ($opts->{minsize} bp)
  --max         integer         Maximum paired-end size allowed ($opts->{maxsize} bp)
  
 Duplication filtering
  --nodedup                     Skip deduplication and use all primary, nondup alignments
  --dupfrac     float           Target duplication rate for subsampling ($opts->{dupfrac})
  --maxdepth    integer         Maximum position alignment depth ($opts->{maxdepth})
                                  set to 1 to remove all duplicates
  --optdist     integer         Maximum distance for optical duplicates ($opts->{optdist})
                                  use 100 for HiSeq, 2500 for NovaSeq
  --deduppair                   Run deduplication as paired-end, but coverage as single-end

 Fragment coverage
  --cutsite                     Set multiple options specific for ATACSeq cutsite analysis
  --size        integer         Predicted fragment size. REQUIRED for single-end
  --shift       integer         Shift the fragment in special situations
  --fraction                    Record multiple-hit alignments as fraction of hits
  --slocal      integer         Small local lambda size ($opts->{slocal} bp)
  --llocal      integer         Large local lambda size ($opts->{llocal} bp)
  --cbin        integer         ChIP fragment bin size ($opts->{chipbin} bp)
  --slbin       integer         Small local lambda bin size ($opts->{slocalbin} bp)
  --llbin       integer         Large local lambda bin size ($opts->{llocalbin} bp)

 Chromosome-specific normalization
  --chrnorm     float           Specific chromosome normalization factor
  --chrapply    "text"          Apply factor to specified chromosomes via regex
 
 Peak calling
  --independent                 Call peaks independently for each replicate and merge
  --cutoff      number          Threshold q-value for calling peaks ($opts->{cutoff}) 
                                  Higher numbers are more significant, -1*log10(q)
  --peaksize    integer         Minimum peak size to call (2 x fragment size)
                                  Required for paired-end alignments.
  --peakgap     integer         Maximum gap between peaks before merging (1 x size)
  --broad                       Also perform broad (gapped) peak calling
  --broadcut    number          Q-value cutoff for linking broad regions ($opts->{broadcut})
  --broadgap    integer         Maximum link size between peaks in broad calls (4 x size bp)
  --nolambda                    Skip lambda control, compare ChIP directly with control
  --minpeakover integer         Minimum number of overlapping replicate peaks to accept
                                  in final when merging (default n-1, minimum 2)
  --samedepth                   Use same target depth when calculating per sample
                                  q-value enrichment for replicate-mean peaks
  
 Peak scoring
  --binsize     integer         Size of bins in 25 flanking peak bins for profile ($opts->{binsize})
  --targetdepth text            Set method for sequence depth scaling for all count data:
                                  median (default), mean, min
  --rawcounts                   Use unscaled raw counts for re-scoring peaks
  --noplot                      Do not plot figures of results
  
 Job control
  --cpu         integer         Number of CPUs to use per job ($opts->{cpu})
  --job         integer         Number of simultaneous jobs ($opts->{job})
  --dryrun                      Just print the commands without execution
  --noorganize                  Do not organize files into subfolders when finished
  --savebam                     Save de-duplicated bam files
  --savebdg                     Save text bedGraph files

 Application Paths
  --bam2wig     path            ($opts->{bam2wig})
  --bamdedup    path            ($opts->{bamdedup})
  --bedtools    path            ($opts->{bedtools})
  --bw2bdg      path            ($opts->{bw2bdg})
  --data2wig    path            ($opts->{data2wig})
  --getdata     path            ($opts->{getdata})
  --getrel      path            ($opts->{getrel})
  --geteff      path            ($opts->{geteff})
  --intersect   path            ($opts->{intersect})
  --macs        path            ($opts->{macs})
  --manwig      path            ($opts->{manwig})
  --meanbdg     path            ($opts->{meanbdg})
  --peak2bed    path            ($opts->{peak2bed})
  --updatepeak  path            ($opts->{updatepeak})
  --pandoc      path            ($opts->{pandoc})
  --plotpeak    path            ($opts->{plotpeak})
  --printchr    path            ($opts->{printchr})
  --reportmap   path            ($opts->{reportmap})
  --rscript     path            ($opts->{rscript})
  --samtools    path            ($opts->{samtools})
  --wig2bw      path            ($opts->{wig2bw})
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
my $argument_string = join q( ), @ARGV;    # save to print below

GetOptions(
	$opts,
	'dir=s',
	'out=s',
	'chip=s@',
	'control=s@',
	'name=s@',
	'chscale=s@',
	'coscale=s@',
	'genome=i',
	'mapq=i',
	'paired|pe!',
	'fraction!',
	'minsize|min=i',
	'maxsize|max=i',
	'chrskip=s',
	'exclude|blacklist=s',
	'dedup!',
	'dupfrac=f',
	'maxdepth=i',
	'optdist=i',
	'deduppair!',
	'fragsize|size=i',
	'shiftsize|shift=i',
	'cutsite|atac!',
	'slocal=i',
	'llocal=i',
	'chipbin|cbin=i',
	'slocalbin|slbin=i',
	'llocalbin|llbin=i',
	'chrnorm=f@',
	'chrapply=s',
	'cutoff=f',
	'targetdepth|tdep=s',
	'samedepth!',
	'peaksize=i',
	'peakgap=i',
	'broad!',
	'broadcut=f',
	'broadgap=i',
	'lambda!',
	'independent!',
	'minpeakover=i',
	'binsize=i',
	'rawcounts!',
	'plot!',
	'cpu=i',
	'job=i',
	'dryrun!',
	'organize!',
	'plot_log2=f',
	'plot_frag=f',
	'plot_qval=f',
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
	'updatepeak=s',
	'combrep=s',
	'plotpeak=s',
	'rscript=s',
	'reportmap=s',
	'samtools=s',
	'pandoc=s'
) or die "unrecognized option(s)!\n";

if ( $opts->{help} ) {
	print $documentation;
	exit 1;
}

### Begin main pipeline
print "======== ChIPSeq multi-replicate pipeline ==========\n";
print "\nversion $VERSION\n";
my $start = time;
check_options();
print_start();
check_input_files();
add_jobs_to_runner();

# Check progress
$Runner->check_progress_file();

# Start pipeline
$Runner->run_generate_chr_file();
$Runner->run_mappable_space_report();
$Runner->run_input_peak_detection();
$Runner->run_dedup();
$Runner->run_bam_filter();
$Runner->run_bam_check();
$Runner->run_bam_fragment_conversion();
$Runner->run_bam_count_conversion();
$Runner->run_bw_conversion();
$Runner->run_lambda_control();
$Runner->run_bdgcmp();
$Runner->run_call_peaks();
$Runner->run_peak_merge();
$Runner->run_bdg_conversion();
$Runner->run_update_peaks();
$Runner->write_samples_file();
$Runner->run_rescore();
$Runner->run_efficiency();
$Runner->run_mean_merge_compare();
$Runner->run_plot_peaks();
$Runner->run_cleanup($argument_string);
$Runner->run_organize();
$Runner->generate_report($argument_string);

# final statement
printf "\n\nFinished in %.1f minutes\n", ( time - $start ) / 60;
print "======== Finished ChIPSeq multi-replicate pipeline ==========\n";

############### Subroutines ########################################################

sub check_options {

	# check chip samples
	if (@ARGV) {
		printf(
"There are unrecognized leftover items on the command line!\n Did you leave spaces in your --chip or --control file lists?\nItems:\n %s\n",
			join( "\n ", @ARGV ) );
		exit 1;
	}
	unless ( $Runner->chips ) {
		print "No ChIP file(s) defined!\n";
		exit 1;
	}
	if ( scalar( $Runner->chips ) == 1 ) {
		print <<END;
This pipeline is intended to work with multiple ChIP samples, but only one was provided.
If you have replicates, try setting each replicate as a separate sample.
Otherwise, try running MACS2 by itself.
END
		exit 1;
	}

	# check names
	if ( $Runner->names ) {
		my %check = map { $_ => 1 } ( $Runner->names );
		if ( scalar( keys %check ) != scalar( $Runner->names ) ) {
			print "Duplicate sample names are present!\n";
			exit 1;
		}
	}
	else {
		print "No name(s) defined!\n";
		exit 1;
	}
	unless ( scalar( $Runner->chips ) == scalar( $Runner->names ) ) {
		print "Unequal number of ChIP samples and names!\n";
		exit 1;
	}

	# check controls
	if (    scalar( $Runner->controls ) > 1
		and scalar( $Runner->controls ) != scalar( $Runner->chips ) )
	{
		print <<END;
Unequal number of control and ChIP samples! Are you missing one?
You may duplicate control file names and the pipeline will handle these smartly.
END
		exit 1;
	}
	elsif ( scalar( $Runner->controls ) == 0 ) {

		# no controls, turn off lambda
		$Runner->slocal(0);
		$Runner->llocal(0);
		$Runner->lambda(0);
	}
	if ( $Runner->slocal == 0 and $Runner->llocal == 0 ) {

		# user somehow set both to zero, but this throws errors to macs2
		$Runner->lambda(0);
	}

	# check scaling factors
	if ( scalar( $Runner->chscales ) or scalar( $Runner->coscales ) ) {

		# no longer recommended
		print <<MESSAGE;

WARNING: Manually setting ChIP and/or Control scaling factors!!!!
Please be aware that manually scaling coverage depth is an advanced option
and should only be done when you're aware of the ramifications. Inappropriate
scaling may artificially alter the expected statistics required for accurate
peak calling and may reduce confidence of identified peaks.

MESSAGE
	}
	if (    scalar( $Runner->chscales )
		and scalar( $Runner->chscales ) != scalar( $Runner->chips ) )
	{
		print "Unequal number of ChIP samples and ChIP scale factors!\n";
		exit 1;
	}
	if (    scalar( $Runner->coscales )
		and scalar( $Runner->coscales ) != scalar( $Runner->controls ) )
	{
		print "Unequal number of control samples and control scale factors!\n";
		exit 1;
	}

	# check chromosome normalization factors
	if ( scalar( $Runner->chrnorm ) and not $Runner->chrapply ) {
		print
"Must specify chromosome apply regex (--chrapply) for chromosome normalization factors!\n";
		exit 1;
	}
	if (    scalar( $Runner->chrnorm )
		and scalar( $Runner->chrnorm ) != scalar( $Runner->chips ) )
	{
		# apply to all the ChIPs
		if ( scalar( $Runner->chrnorm ) == 1 ) {
			print "Using the same chromosome normalization factor for each ChIP sample\n";
			my $n = ( $Runner->chrnorm )[0];
			my $m = scalar( $Runner->names );
			for ( 2 .. $m ) {

				# add it to the remainder
				$Runner->chrnorm($n);
			}
		}
		else {
			print "Unequal number of chromosome normalization factors and ChIP samples!\n";
			exit 1;
		}
	}
	if ( scalar( $Runner->chrnorm ) > 1 and scalar( $Runner->controls ) == 1 ) {
		print "Using first chromosome normalization factor for universal control!\n";
	}

	# check directory
	unless ( $Runner->dryrun ) {
		unless ( -e $Runner->dir and -d _ ) {
			make_path( $Runner->dir )
				or die
				sprintf( "Unable to make directory %s! $OS_ERROR\n", $Runner->dir );
		}
		unless ( -w $Runner->dir ) {
			die sprintf( "Target directory %s is not writable!\n", $Runner->dir );
		}
	}

	# check target depth
	if ( $Runner->targetdepth and $Runner->targetdepth =~ /^ \d+ (?: \.\d+)? $/x) {
		my $files = join( ",", ( $Runner->chips ), ( $Runner->controls ) );
		if ( $files =~ /\.bam/i ) {

			# no longer recommend manually setting with bam files
			# bigWig files would be acceptable
			print <<MESSAGE;

WARNING: Manually setting the target sequence depth is not recommended!!!!
Target depth is now automatically calculated empirically (see documentation).
Reconsider using this option unless you understand the ramifications.

MESSAGE
		}
	}
	elsif ( $Runner->targetdepth ) {
		if ( $Runner->targetdepth !~ /mean | median | min/x ) {
			printf("unrecognized targetdepth parameter '%s'!", $Runner->targetdepth);
			exit 1;
		}
	}
	else {
		$Runner->targetdepth('median');
	}

	# check sizes
	if ( $Runner->cutsite ) {
		# set special options for ATACSeq cutsite analysis
		print 
" Using specific parameters for ATAC Cut Site mode\n These may override manual settings\n";
		$Runner->deduppair(1);
		$Runner->paired(0);
		$Runner->minsize(30);    # minsize and maxsize technically not required
		$Runner->maxsize(2000);
		$Runner->fragsize(50);
		$Runner->shiftsize(-25);
		$Runner->peaksize(90);
		$Runner->peakgap(30);
		$Runner->slocal(100);    # very rare to use lambda with ATAC but just in case
		$Runner->llocal(1000);
		$Runner->broad(0);
	}
	if ( not $Runner->peaksize ) {

		# no minimum peak size defined? might be ok
		if ( $Runner->fragsize ) {

			# set default to twice fragment size
			$Runner->peaksize( 2 * $Runner->fragsize );
		}
		elsif ( $Runner->paired ) {

			# paired fragments - make something up
			if ( not $Runner->peaksize ) {
				print
"\nWARNING! Setting minimum peak size to 500 bp, but this should be manually\nset based on mean alignment insert size and nature of experiment.\n\n";
				$Runner->peaksize(500);
				$Runner->fragsize(250);    # for expected background normalization
			}
		}
		else {
			print <<END;
 Must set an estimated mean fragment size for single-end alignments!
 Run 'macs2 predictd' or 'bam2wig.pl --shift --model'
END
			exit 1;
		}
	}
	unless ( defined $Runner->peakgap ) {
		my $g = $Runner->fragsize ? $Runner->fragsize : int( $Runner->peaksize / 2 );
		$Runner->peakgap($g);
	}
	unless ( defined $Runner->broadgap ) {
		my $g = $Runner->fragsize ? 4 * $Runner->fragsize : 2 * $Runner->peaksize;
		$Runner->broadgap($g);
	}

	# check exclusion list
	if ( my $exclude = $Runner->exclude ) {
		if ( $exclude eq 'none' ) {
			# this is ok, we ignore
		}
		elsif ( $exclude eq 'input' ) {
			# this is ok so long as we actually have controls
			unless ( scalar( $Runner->controls ) ) {
				print
				"\nWARNING! No Control files, cannot use 'input' as exclusion list!\n";
				$Runner->exclude('none');
			}
		}
		else {
			unless ( -e $exclude ) {
				print
"\nWARNING! Unable to find specified exclusion list file '%s'!\n Defaulting to 'none'\n",
				$exclude;
				$Runner->exclude('none');
			}
		}
	}
	elsif ( not defined $Runner->exclude and scalar( $Runner->controls ) ) {
		$Runner->exclude('input');
	}

	# max depth-duplication confusion
	if ( $Runner->maxdup ) {

		# because this was inappropriately named before
		print " \nWARNING: The --maxdup option is now --maxdepth\n";
		$Runner->maxdepth( $Runner->maxdup );
	}
}

sub print_start {
	if ( $Runner->dryrun ) {
		print <<DRYRUN;

======= Dry Run =======
Some values are empirically determined during execution and are made up here.
Some files may not be generated during actual execution, but commands should
mostly be complete. No files will be generated.
=======================

DRYRUN
	}
	print "\nProvided options: $argument_string\n\n";
	$Runner->print_config;
	print "\n\n";
}

sub check_input_files {
	print "\n\n======= Checking input files\n";
	my @errors;
	my @names          = $Runner->names;
	my @chips          = $Runner->chips;
	my @controls       = $Runner->controls;
	my $chip_max_count = 0;
	foreach my $i ( 0 .. $#names ) {
		my @list = split /,/, $chips[$i];
		if ( scalar(@list) != uniqstr(@list) ) {
			printf " Fixing duplicate entries found in %s ChIP entries!\n", $names[$i];

			# fix it
			@list = uniqstr(@list);

			# hack to replace
			$chips[$i] = join( ',', @list );
			$opts->{chip} = \@chips;
		}
		foreach my $f (@list) {
			unless ( -e $f ) {
				printf( " can't find %s ChIP file %s!\n", $names[$i], $f );
				push @errors, $f;
			}
		}
		if ( scalar(@list) > $chip_max_count ) {
			$chip_max_count = scalar(@list);
		}

		if (@controls) {
			my @list2 = split /,/, $controls[$i];
			if ( scalar(@list2) != uniqstr(@list2) ) {
				printf " Fixing duplicate entries found in %s control entries!\n",
					$names[$i];

				# fix it
				@list2 = uniqstr(@list2);

				# hack to replace
				$controls[$i] = join( ',', @list2 );
				$opts->{control} = \@controls;
			}
			foreach my $f (@list2) {
				unless ( -e $f ) {
					printf( " can't find %s control file %s!\n", $names[$i], $f );
					push @errors, $f;
				}
			}
		}
	}
	if ( scalar @errors ) {
		die sprintf(
			"ERROR! Missing %d input files!\nMissing: %s\n",
			scalar @errors, join( ', ', @errors )
		);
	}
	else {
		print " All input files found\n";
	}
	if ( $Runner->independent and $chip_max_count == 1 ) {
		print " Disabling independent flag because no sample has more than 1 file\n";
		$Runner->independent(0);
	}
}

sub add_jobs_to_runner {

	# we pass this off to ChIPjob package with 6 options:
	# name, chip files, control files, chip scale factor, control scale factor,
	# chromosome normalization factor

	# first check for a unversal control
	if ( scalar( $Runner->controls ) == 1 and scalar( $Runner->chips ) > 1 ) {
		print " Using only one control for multiple ChIP experiments\n";
		my $universal_control = ( $Runner->controls )[0];
		$Runner->has_universal_control(1);

		# generate universal name
		my $universal_name;
		if ( $universal_control =~ /\.bam$/i ) {

			# one bam file
			$universal_name = $Runner->out . '_control';
			if ( $Runner->lambda ) {

				# add lambda control if we're using that, otherwise leave it be
				$universal_name .= '.lambda_control';
			}
			else {
				$universal_name .= '.control_fragment';
			}
		}
		elsif ( $universal_control =~ /\. (?: bw | bigwig )$/xi ) {

			# a pre-processed bigWig file
			( undef, undef, $universal_name ) = File::Spec->splitpath($universal_control);
			$universal_name =~ s/\. (?: bw | bigwig )$//x;
		}
		else {
			die "unrecognized control file '$universal_control'!\n";
		}
		my $universal_scale = ( $Runner->coscales )[0] || undef;

		# add back universal control with special prefix for each ChIP job
		$opts->{control} = [];    # unorthodox hack
		$opts->{coscale} = [];    # unorthodox hack
		foreach ( $Runner->names ) {
			$Runner->controls( 'Custom-Universal-' . $universal_name );
		}

		# generate job - this will always be the first job
		$Runner->add_job(
			$universal_name,
			q(),
			$universal_control,
			undef,
			$universal_scale,
			( $Runner->chrnorm )[0] || undef
		);
	}

	# walk through each given job
	my @names = $Runner->names;
	for my $i ( 0 .. $#names ) {

		# add job
		$Runner->add_job(
			$names[$i],
			( $Runner->chips )[$i],
			( $Runner->controls )[$i] || undef,
			( $Runner->chscales )[$i] || undef,
			( $Runner->coscales )[$i] || undef,
			( $Runner->chrnorm )[$i] || undef
		);
	}
}

