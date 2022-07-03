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
use File::Spec;
use File::Path qw(make_path);
use Getopt::Long;
use Bio::MultiRepChIPSeq::Runner;


# Initialize Runner and options
my $Runner = Bio::MultiRepChIPSeq::Runner->new();
my $opts = $Runner->options;
my $VERSION = $Runner->version;

# reset some default parameters for recalling
undef $opts->{peaksize};
undef $opts->{peakgap};
undef $opts->{broadcut};
undef $opts->{broadgap};

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

New peak and analysis files are placed into new subdirectories with a 
numeric suffix (allowing for multiple conditions to be run consecutively) 
without overwriting pre-existing files. Peaks from different runs can 
subsequently be compared using the intersect_peaks.pl script, if desired.

See https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq for full 
usage and guide.

Version: $VERSION

Options:
 Input files
  --in        file basename     Base filename for previous pipeline output
  --dir       directory         Directory for writing all files ($opts->{dir})
  --out       file basename     Base filename for new output files ($opts->{out})
  
 Peak calling
  --cutoff    number            Threshold q-value for calling peaks () 
                                  Higher numbers are more significant, -1*log10(q)
  --peaksize  integer           Minimum peak size to call ()
  --peakgap   integer           Maximum gap between peaks before merging ()
  --broad                       Also perform broad (gapped) peak calling
  --broadcut  number            Q-value cutoff for linking broad regions ()
  --broadgap  integer           Maximum link size between peaks in broad calls ()
  
 Peak scoring
  --binsize   integer           Size of bins in 25 flanking peak bins for profile ($opts->{binsize})
  --repmean                     Combine replicate counts as mean for each sample set
  --noplot                      Do not plot figures of results
  
 Job control
  --cpu       integer           Number of CPUs to use per job ($opts->{cpu})
  --job       integer           Number of simultaneous jobs ($opts->{job})
  --dryrun                      Just print the commands without execution
  --noorganize                  Do not organize files into subfolders when finished

 Application  Paths
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
my $argument_string = join ' ', @ARGV; # save to print below

GetOptions(
	$opts,
	'in=s',
	'dir=s',
	'out=s',
	'cutoff=f',
	'peaksize=i',
	'peakgap=i',
	'broad!',
	'broadcut=f',
	'broadgap=i',
	'plot!',
	'cpu=i',
	'job=i',
	'dryrun!',
	'organize!',
	'help!',
	'chromofile=s',
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
$Runner->run_clean_peaks();
$Runner->run_peak_merge();
$Runner->run_bdg_conversion();
$Runner->run_rescore();
$Runner->run_efficiency();
$Runner->run_plot_peaks();
$Runner->run_cleanup($argument_string);
$Runner->run_organize($suffix);

# final statement
printf "\n\nFinished in %.1f minutes\n", (time -$start) / 60;
print "======== Finished ChIPSeq Peak Re-Caller ==========\n";




############### Subroutines ########################################################

sub check_inputs {
	if (@ARGV) {
		die sprintf("There are unrecognized leftover items on the command line!\n Did you leave spaces in your --chip or --control file lists?\nItems:\n %s\n", join("\n ", @ARGV));
	}
	
	# sizes
	unless ($Runner->cutoff and $Runner->peaksize and $Runner->peakgap) {
		die "Must set --cutoff, --peaksize, and --peakgap paramets!\n";
	}
	$Runner->broad(1) if ($Runner->broadgap or $Runner->broadcut);
	if ($Runner->broad) {
		die "Must set --broadcut if doing broad (gapped) peak calling!\n" 
			unless $Runner->broadcut;
		die "Must set --broadgap if doing broad (gapped) peak calling!\n" 
			unless $Runner->broadgap;
	}
	
	# input / output
	unless ($Runner->in) {
		die "Must set the previous pipeline run output basename using --in\n";
	}
	unless ($Runner->out) {
		die "Must set a new output basename with --out\n";
	}
	if ($Runner->in eq $Runner->out) {
		die "Please set different input and output base names\n";
	}
	
	# check suffix
	if ($Runner->organize) {
		my $i = 1;
		my $check = File::Spec->catfile($Runner->dir, 'Peaks' . $i);
		while (-e $check) {
			$i++;
			$check = File::Spec->catfile($Runner->dir, 'Peaks' . $i);
		}
		$suffix = $i;
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
	print "\nProvided options: $argument_string\n\n";
	$Runner->print_config;
	print "\n\n";
}

sub check_input_files {
	print "\n\n======= Checking input files\n";
	
	# find sample file list
	my $sample_file = File::Spec->catfile($Runner->dir, $Runner->in . '_samples.txt');
	my $organized = 0;
	unless (-e $sample_file) {
		# it may have been organized
		my $sf2 = File::Spec->catfile($Runner->dir, 'Analysis', $Runner->in . '_samples.txt');
		if (-e $sf2) {
			$sample_file = $sf2;
			$organized = 1;
		}
		else {
			die sprintf("Unable to identify sample file '%s' in either directories %s or %s!\n", 
				$Runner->in . '_samples.txt', $Runner->dir, 
				File::Spec->catfile($Runner->dir, 'Analysis'));
		}
	}
	
	# load samples
	my $fh = IO::File->new($sample_file) or 
		die "unable to read $sample_file! $!\n";
	my $h = $fh->getline;
	unless ($h eq "Replicate\tDataset\n") {
		die "Sample file '$sample_file' does not have the requisite header!\n";
	}
	my %samples;
	my @sample_order;
	while (my $line = $fh->getline) {
		chomp $line;
		my ($replicate, $condition) = split('\t', $line);
		unless (exists $samples{$condition}) {
			$samples{$condition} = [];
			if (lc $condition ne 'input') {
				push @sample_order, $condition;
			}
		}
		push @{ $samples{$condition} }, $replicate;
	}
	$fh->close;
	
	# Add jobs and associated files for each condition
	my @errors;
	foreach my $condition (@sample_order) {
		my $Job = $Runner->add_job($condition, undef, undef, undef, undef, undef);
			# Runner will crash if Job wasn't made succesfully
			
		# Add existing track files 
		if ($organized) {
			# files based on known organizational subfolders
			
			my $qbw = File::Spec->catfile($Runner->dir, 'QValue', $condition . '.qvalue.bw');
			if (-r $qbw) {
				$Job->qvalue_bw($qbw);
				$Job->qvalue_bdg( 
					File::Spec->catfile($Runner->dir, $condition . '.qvalue.bdg')
				);
			}
			else {
				push @errors, "  File $qbw cannot be found!\n";
			}
			
			my $lbw = File::Spec->catfile($Runner->dir, 'Log2FE', $condition . '.log2FE.bw');
			if (-r $lbw) {
				$Job->logfe_bw($lbw);
			}
			else {
				push @errors, "  File $lbw cannot be found!\n";
			}
			
			my $fbw = File::Spec->catfile($Runner->dir, 'Fragment', $condition . '.fragment.bw');
			if (-r $fbw) {
				$Job->chip_bw($fbw);
			}
			else {
				push @errors, "  File $fbw cannot be found!\n";
			}
			
			my @cbw = map { File::Spec->catfile($Runner->dir, 'Count', $_ . '.count.bw') }
				@{ $samples{$condition} };
			foreach my $c (@cbw) {
				if (-r $c) {
					$Job->chip_count_bw($c);
				}
				else {
					push @errors, "  File $c cannot be found!\n";
				}
			}
			
			if (exists $samples{'Input'}) {
				# I won't necessarily know which is which, so do all of them????
				my @icb = map {
					File::Spec->catfile($Runner->dir, 'Count', $_ . '.count.bw') 
					} @{ $samples{'Input'} };
				foreach my $c (@icb) {
					if (-r $c) {
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

			my $qbw = File::Spec->catfile($Runner->dir, $condition . '.qvalue.bw');
			if (-r $qbw) {
				$Job->qvalue_bw($qbw);
				$Job->qvalue_bdg(
					File::Spec->catfile($Runner->dir, $condition . '.qvalue.bdg')
				);
			}
			else {
				push @errors, "  File $qbw cannot be found!\n";
			}
			
			my $lbw = File::Spec->catfile($Runner->dir, $condition . '.log2FE.bw');
			if (-r $lbw) {
				$Job->logfe_bw($lbw);
			}
			else {
				push @errors, "  File $lbw cannot be found!\n";
			}
			
			my $fbw = File::Spec->catfile($Runner->dir, $condition . '.fragment.bw');
			if (-r $fbw) {
				$Job->chip_bw($fbw);
			}
			else {
				push @errors, "  File $fbw cannot be found!\n";
			}
			
			my @cbw = map { File::Spec->catfile($Runner->dir, $_ . '.count.bw') }
				@{ $samples{$condition} };
			foreach my $c (@cbw) {
				if (-r $c) {
					$Job->chip_count_bw($c);
				}
				else {
					push @errors, "  File $c cannot be found!\n";
				}
			}
			
			if (exists $samples{'Input'}) {
				# I won't necessarily know which is which, so do all of them????
				my @icb = map { File::Spec->catfile($Runner->dir, $_ . '.count.bw') }
					@{ $samples{'Input'} };
				foreach my $c (@icb) {
					if (-r $c) {
						$Job->control_count_bw($c);
					}
					else {
						push @errors, "  File $c cannot be found!\n";
					}
				}
			}
		}
		
		# add new peak files
		$Job->peak(File::Spec->catfile($Runner->dir, $condition . '.narrowPeak'));
		$Job->clean_peak(File::Spec->catfile($Runner->dir, $condition . '.bed'));
		$Job->peak_summit(File::Spec->catfile($Runner->dir, $condition . '.summit.bed'));
		if ($Runner->broad) {
			$Job->gappeak(File::Spec->catfile($Runner->dir, $condition . '.gappedPeak'));
			$Job->clean_gappeak(File::Spec->catfile($Runner->dir, $condition . '.gapped.bed'));
		}
		
		# unset presumed lambda files
		$Job->lambda_bw(q());
		$Job->lambda_bdg(q());
	}
	
	# summary
	foreach my $Job ($Runner->list_jobs) {
		printf " ChIP Job '%s' with %d ChIP replicates and %d Input replicates\n",
			$Job->job_name, scalar($Job->chip_count_bw), scalar($Job->control_count_bw);
	}
	print "\n\n";
}



