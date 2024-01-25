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
use IO::File;
use File::Spec;
use File::Path qw(make_path);
use Getopt::Long;
use Bio::MultiRepChIPSeq::Runner;

# Initialize Runner and options
my $Runner = Bio::MultiRepChIPSeq::Runner->new();
my $opts   = $Runner->options;
our $VERSION = $Runner->version;

# reset some default parameters for recalling
undef $opts->{peaksize};
undef $opts->{peakgap};
undef $opts->{broadcut};
undef $opts->{broadgap};
undef $opts->{dedup};

my $documentation = <<DOC;

======= ChIPSeq Peak Re-Caller ==========

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

NOTE: Independently called peaks for replicates are not re-called. 
Only mean-replicate peaks are recalled. 

New peak and analysis files are placed into new subdirectories with a 
numeric suffix (allowing for multiple conditions to be run consecutively) 
without overwriting pre-existing files. Peaks from different runs can 
subsequently be compared using the intersect_peaks.pl script, if desired.

See https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq for full 
usage and guide.

Version: $VERSION

Options:
 Input files
  --in          file basename   Base filename for previous pipeline output
  --dir         directory       Directory for writing all files ($opts->{dir})
  --out         file basename   Base filename for new output files ($opts->{out})
  
 Peak calling
  --cutoff      number          Threshold q-value for calling peaks () 
                                  Higher numbers are more significant, -1*log10(q)
  --peaksize    integer         Minimum peak size to call ()
  --peakgap     integer         Maximum gap between peaks before merging ()
  --broad                       Also perform broad (gapped) peak calling
  --broadcut    number          Q-value cutoff for linking broad regions ()
  --broadgap    integer         Maximum link size between peaks in broad calls ()
  --atac                        Convenience option to set peak size and gap
  
 Peak scoring
  --binsize     integer         Size of bins in 25 flanking peak bins for profile ($opts->{binsize})
  --noplot                      Do not plot figures of results
  
 Job control
  --cpu         integer         Number of CPUs to use per job ($opts->{cpu})
  --job         integer         Number of simultaneous jobs ($opts->{job})
  --dryrun                      Just print the commands without execution
  --noorganize                  Do not organize files into subfolders when finished

 Application Paths
  --macs        path            ($opts->{macs})
  --manwig      path            ($opts->{manwig})
  --wig2bw      path            ($opts->{wig2bw})
  --bw2bdg      path            ($opts->{bw2bdg})
  --printchr    path            ($opts->{printchr})
  --data2wig    path            ($opts->{data2wig})
  --getdata     path            ($opts->{getdata})
  --getrel      path            ($opts->{getrel})
  --geteff      path            ($opts->{geteff})
  --meanbdg     path            ($opts->{meanbdg})
  --bedtools    path            ($opts->{bedtools})
  --intersect   path            ($opts->{intersect})
  --pandoc      path            ($opts->{pandoc})
  --peak2bed    path            ($opts->{peak2bed})
  --updatepeak  path            ($opts->{updatepeak})
  --plotpeak    path            ($opts->{plotpeak})
  --rscript     path            ($opts->{rscript})
DOC

### Inputs
unless (@ARGV) {
	print <<HELP;

======= ChIPSeq Peak Re-Caller ==========

This is a script for re-calling the peaks with new parameters from the 
output of a the MultiRep ChIPSeq Pipeline. 

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
	'in=s',
	'dir=s',
	'out=s',
	'atac!',
	'cutoff=f',
	'peaksize=i',
	'peakgap=i',
	'broad!',
	'broadcut=f',
	'broadgap=i',
	'binsize=i',
	'plot!',
	'cpu=i',
	'job=i',
	'dryrun!',
	'organize!',
	'plot_log2=f',
	'plot_frag=f',
	'plot_qval=f',
	'help!',
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
	'plotpeak=s',
	'rscript=s',
	'pandoc=s'
) or die "unrecognized option(s)!\n";

if ( $opts->{help} ) {
	print $documentation;
	exit 1;
}

### Begin main pipeline
print "======== ChIPSeq Peak Re-Caller ==========\n";
print "\nversion $VERSION\n";
my $start = time;
my $suffix;
check_inputs();
print_start();
check_input_files();

# Start pipeline
$Runner->check_progress_file();
$Runner->run_generate_chr_file();
$Runner->run_bw_conversion();
$Runner->run_call_peaks();
$Runner->run_peak_merge();
$Runner->run_bdg_conversion();
$Runner->run_update_peaks();
$Runner->run_rescore();
$Runner->run_efficiency();
$Runner->run_plot_peaks();
$Runner->run_cleanup($argument_string);
$Runner->run_organize($suffix);
$Runner->generate_report($argument_string);

# final statement
printf "\n\nFinished in %.1f minutes\n", ( time - $start ) / 60;
print "======== Finished ChIPSeq Peak Re-Caller ==========\n";

############### Subroutines ########################################################

sub check_inputs {
	if (@ARGV) {
		printf(
"There are unrecognized leftover items on the command line!\n Did you leave spaces in your --chip or --control file lists?\nItems:\n %s\n",
			join( "\n ", @ARGV ) );
		exit 1;
	}

	# sizes
	if ( $Runner->atac ) {
		# just to maintain consistency with main pipeline
		# only useful parameters is to set peak size and gap
		unless ( $Runner->peaksize ) {
			$Runner->peaksize(150);
		}
		unless ( $Runner->peakgap ) {
			$Runner->peakgap(50);
		}
	}
	unless ( $Runner->cutoff and $Runner->peaksize and $Runner->peakgap ) {
		print "Must set --cutoff, --peaksize, and --peakgap paramets!\n";
		exit 1;
	}
	$Runner->broad(1) if ( $Runner->broadgap or $Runner->broadcut );
	if ( $Runner->broad ) {
		unless ( $Runner->broadcut ) {
			print "Must set --broadcut if doing broad (gapped) peak calling!\n";
			exit 1;
		}
		unless ( $Runner->broadgap ) {
			print "Must set --broadgap if doing broad (gapped) peak calling!\n";
			exit 1;
		}
	}

	# input / output
	unless ( $Runner->in ) {
		print "Must set the previous pipeline run output basename using --in\n";
		exit 1;
	}
	unless ( $Runner->out ) {
		print "Must set a new output basename with --out\n";
		exit 1;
	}
	if ( $Runner->in eq $Runner->out ) {
		print "Please set different input and output base names\n";
		exit 1;
	}

	# check suffix
	if ( $Runner->organize ) {
		my $i     = 1;
		my $check = File::Spec->catfile( $Runner->dir, 'Peaks' . $i );
		while ( -e $check ) {
			$i++;
			$check = File::Spec->catfile( $Runner->dir, 'Peaks' . $i );
		}
		$Runner->dir_suffix($i);
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

	# find sample file list
	my $dir       = File::Spec->catfile( $Runner->dir, 'Analysis' );
	my $organized = 0;
	if ( -e $dir ) {
		$organized = 1;
	}
	else {
		$dir = $Runner->dir;
		if ( -e $dir ) {
			$organized = 0;
		}
		else {
			die " Directory $dir does not exist!?";
		}
	}
	my $sample_file = File::Spec->catfile( $dir, $Runner->in . '_samples.txt' );
	unless ( -e $sample_file ) {
		$sample_file = File::Spec->catfile( $dir, $Runner->in . '.rep_mean_samples.txt' );
		if ( not -e $sample_file ) {
			die " Unable to identify sample files in $dir!\n";
		}
	}

	# load samples
	my $fh = IO::File->new($sample_file)
		or die "unable to read $sample_file! $OS_ERROR\n";
	my $h = $fh->getline;
	unless ( $h eq "Replicate\tDataset\n" ) {
		die "Sample file '$sample_file' does not have the requisite header!\n";
	}
	my %samples;
	my @sample_order;
	while ( my $line = $fh->getline ) {
		chomp $line;
		my ( $replicate, $condition ) = split /\t/, $line;
		unless ( exists $samples{$condition} ) {
			$samples{$condition} = [];
			if ( lc $condition ne 'input' ) {
				push @sample_order, $condition;
			}
		}
		push @{ $samples{$condition} }, $replicate;
	}
	$fh->close;

	# Add jobs and associated files for each condition
	my @errors;
	foreach my $condition (@sample_order) {
		my $Job = $Runner->add_job( $condition, undef, undef, undef, undef, undef );

		# Runner will crash if Job wasn't made succesfully

		# Add existing track files
		if ($organized) {

			# files based on known organizational subfolders

			my $qbw =
				File::Spec->catfile( $Runner->dir, 'QValue', $condition . '.qvalue.bw' );
			if ( -r $qbw ) {
				$Job->qvalue_bw($qbw);
				$Job->qvalue_bdg(
					File::Spec->catfile( $Runner->dir, $condition . '.qvalue.bdg' ) );
			}
			else {
				push @errors, "  File $qbw cannot be found!\n";
			}

			my $lbw =
				File::Spec->catfile( $Runner->dir, 'Log2FE', $condition . '.log2FE.bw' );
			if ( -r $lbw ) {
				$Job->logfe_bw($lbw);
			}
			else {
				push @errors, "  File $lbw cannot be found!\n";
			}

			my $fbw = File::Spec->catfile( $Runner->dir, 'Fragment',
				$condition . '.fragment.bw' );
			if ( -r $fbw ) {
				$Job->chip_bw($fbw);
			}
			else {
				push @errors, "  File $fbw cannot be found!\n";
			}

			my @cbw =
				map { File::Spec->catfile( $Runner->dir, 'Count', $_ . '.count.bw' ) }
				@{ $samples{$condition} };
			foreach my $c (@cbw) {
				if ( -r $c ) {
					$Job->chip_count_bw($c);
				}
				else {
					push @errors, "  File $c cannot be found!\n";
				}
			}

			if ( exists $samples{'Input'} ) {

				# I won't necessarily know which is which, so do all of them????
				my @icb =
					map { File::Spec->catfile( $Runner->dir, 'Count', $_ . '.count.bw' ) }
					@{ $samples{'Input'} };
				foreach my $c (@icb) {
					if ( -r $c ) {
						$Job->control_count_bw($c);
					}
					else {
						push @errors, "  File $c cannot be found!\n";
					}
				}
			}
		}
		else {
			# no known organizational subfolders

			my $qbw = File::Spec->catfile( $Runner->dir, $condition . '.qvalue.bw' );
			if ( -r $qbw ) {
				$Job->qvalue_bw($qbw);
				$Job->qvalue_bdg(
					File::Spec->catfile( $Runner->dir, $condition . '.qvalue.bdg' ) );
			}
			else {
				push @errors, "  File $qbw cannot be found!\n";
			}

			my $lbw = File::Spec->catfile( $Runner->dir, $condition . '.log2FE.bw' );
			if ( -r $lbw ) {
				$Job->logfe_bw($lbw);
			}
			else {
				push @errors, "  File $lbw cannot be found!\n";
			}

			my $fbw = File::Spec->catfile( $Runner->dir, $condition . '.fragment.bw' );
			if ( -r $fbw ) {
				$Job->chip_bw($fbw);
			}
			else {
				push @errors, "  File $fbw cannot be found!\n";
			}

			my @cbw = map { File::Spec->catfile( $Runner->dir, $_ . '.count.bw' ) }
				@{ $samples{$condition} };
			foreach my $c (@cbw) {
				if ( -r $c ) {
					$Job->chip_count_bw($c);
				}
				else {
					push @errors, "  File $c cannot be found!\n";
				}
			}

			if ( exists $samples{'Input'} ) {

				# I won't necessarily know which is which, so do all of them????
				my @icb = map { File::Spec->catfile( $Runner->dir, $_ . '.count.bw' ) }
					@{ $samples{'Input'} };
				foreach my $c (@icb) {
					if ( -r $c ) {
						$Job->control_count_bw($c);
					}
					else {
						push @errors, "  File $c cannot be found!\n";
					}
				}
			}
		}

		# add new peak files
		$Job->repmean_peak(
			File::Spec->catfile( $Runner->dir, $condition . '.narrowPeak' ) );
		if ( $Runner->broad ) {
			$Job->repmean_gappeak(
				File::Spec->catfile( $Runner->dir, $condition . '.gappedPeak' ) );
		}

		# unset presumed lambda files
		$Job->lambda_bw(q());
		$Job->lambda_bdg(q());
	}

	# summary
	foreach my $Job ( $Runner->list_jobs ) {
		printf " ChIP Job '%s' with %d ChIP replicates and %d Input replicates\n",
			$Job->job_name, scalar( $Job->chip_count_bw ),
			scalar( $Job->control_count_bw );
	}
	print "\n\n";
}

