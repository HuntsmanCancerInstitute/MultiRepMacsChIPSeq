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
use File::Which;
use File::Path qw(make_path);
use Getopt::Long;

my $VERSION = 12.1;

my $parallel;
eval {
	require Parallel::ForkManager;
	$parallel = 1;
};

# parameters
my %opts = (
	dir         => './',
	out         => 'merged',
	genome      => 0,
	species     => 'human',
	paired      => 0,
	mapq        => 0,
	minsize     => 50,
	maxsize     => 500,
	dup         => 1,
	maxdup      => undef,
	dupfrac     => 0.1,
	optdist     => 0,
	deduppair   => 0,
	savebam     => 0,
	fragsize    => 250,
	shiftsize   => 0,
	slocal      => 1000,
	llocal      => 10000,
	qvalue      => 2,
	peaksize    => undef,
	peakgap     => undef,
	broad       => 0,
	linkqv      => 1,
	gaplink     => undef,
	targetdep   => 25,
	use_lambda  => 1,
	chrskip     => "chrM|MT|lambda|Adapter|PhiX",
	blacklist   => undef,
	cpu         => 4,
	chipbin     => 10,
	slocalbin   => 50,
	llocalbin   => 100,
	chrapply    => undef,
	rawcounts   => 0,
	savebdg     => 0,
	genomewin   => 0,
	discard     => 10,
	repmean     => 0,
	plot        => 0,
	bam2wig     => sprintf("%s", which 'bam2wig.pl'),
	bamdedup    => sprintf("%s", which 'bam_partial_dedup.pl'),
	macs        => sprintf("%s", which 'macs2'),
	manwig      => sprintf("%s", which 'manipulate_wig.pl'),
	wig2bw      => sprintf("%s", which 'wigToBigWig'),
	bw2bdg      => sprintf("%s", which 'bigWigToBedGraph'),
	bedtools    => sprintf("%s", which 'bedtools'),
	getdata     => sprintf("%s", which 'get_datasets.pl'),
	printchr    => sprintf("%s", which 'print_chromosome_lengths.pl'),
	data2wig    => sprintf("%s", which 'data2wig.pl'),
	meanbdg     => sprintf("%s", which 'generate_mean_bedGraph.pl'),
	intersect   => sprintf("%s", which 'intersect_peaks.pl'),
	combrep     => sprintf("%s", which 'combine_replicate_data.pl'),
	plotpeak    => sprintf("%s", which 'plot_peak_figures.R'),
	rscript     => sprintf("%s", which 'Rscript'),
);
$opts{job} = $parallel ? 2 : 1;
my @names;
my @chips;
my @controls;
my @chip_scales;
my @control_scales;
my @chrnorms;

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
If no control file is provided, then the chromosomal mean from the ChIP file is used 
as a (poor) substitute. 

When using a calibration genome in the ChIPSeq (ChIP-Rx), calculate the ratio of 
target/reference alignments for each replica and sample. Provide these with the 
scale options in the same order as the bam files. Note that significant 
de-duplication levels may affect these ratios; if possible, de-duplicate first.

Advanced users may provide one processed bigWig file per ChIP or control sample. 

See https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq for full 
usage and guide.

Version: $VERSION

Options:
 Input files
  --chip      file1,file2...    Repeat for each sample set
  --name      text              Repeat for each sample
  --control   file1,file2...    Repeat if matching multiple samples
  --chscale   number,number...  Calibration scales for each ChIP replica and set
  --coscale   number,number...  Calibration scales for each control replica and set
 
 Output
  --dir       directory         Directory for writing all files ($opts{dir})
  --out       file basename     Base filename for merged output files ($opts{out})
  
 Genome size
  --species   [human,mouse,fish,fly,yeast]   Default ($opts{species})
  --genome    integer           Alternatively give effective genome size
  
 Bam options
  --mapq      integer           Minimum mapping quality, ($opts{mapq})
  --pe                          Bam files are paired-end, default treat as single-end
  --min       integer           Minimum paired-end size allowed ($opts{minsize} bp)
  --max       integer           Maximum paired-end size allowed ($opts{maxsize} bp)
 
 Bam filtering options
  --chrskip   "text"            Chromosome skip regex ($opts{chrskip})
  --blacklist file              Bed file of repeats or hotspots to avoid
  
 Duplication filtering
  --nodup                       Skip deduplication
  --dupfrac   fraction          Minimum allowed fraction of duplicates ($opts{dupfrac})
  --maxdup    integer           Maximum allowed duplication depth ($opts{maxdup})
  --optdist   integer           Maximum distance for optical duplicates ($opts{optdist})
                                  use 100 for HiSeq, 10000 for NovaSeq
  --deduppair                   Run deduplication as paired-end only
  --savebam                     Save de-duplicated bam files

 Fragment coverage
  --size      integer           Predicted fragment size (single-end only, $opts{fragsize} bp)
  --shift     integer           Shift the fragment, e.g. ATACSeq ($opts{shiftsize} bp)
  --slocal    integer           Small local lambda size ($opts{slocal} bp)
  --llocal    integer           Large local lambda size ($opts{llocal} bp)
  --cbin      integer           ChIP fragment bin size ($opts{chipbin} bp)
  --slbin     integer           Small local lambda bin size ($opts{slocalbin} bp)
  --llbin     integer           Large local lambda bin size ($opts{llocalbin} bp)

 Chromosome-specific normalization
  --chrnorm   fraction          Specific chromosome normalization factor
  --chrapply  "text"            Apply factor to specified chromosomes
 
 Peak calling
  --cutoff    number            Threshold q-value for calling peaks ($opts{cutoff}) 
                                 Higher numbers are more significant, -1*log10(q)
  --tdep      integer           Average sequence depth of bam files in millions ($opts{targetdep})
  --peaksize  integer           Minimum peak size to call (2 x size)
  --peakgap   integer           Maximum gap between peaks before merging (1 x size)
  --broad                       Also perform broad (gapped) peak calling
  --broadcut  number            Q-value cutoff for linking broad regions ($opts{broadcut})
  --broadgap  integer           Maximum link size between peaks in broad calls (4 x size bp)
  --nolambda                    Skip lambda control, compare ChIP directly with control
  --rawcounts                   Use unscaled raw counts for re-scoring peaks
  --savebdg                     Save q-value bdg files for further custom calling
  --window    integer           Collect counts across genome in given window size
  --discard   number            Discard genome windows with replicate sum below number ($opts{discard})
  --repmean                     Combine replicate counts as mean for each sample set
  --plot                        Plot figures of results
  
 Job control
  --cpu       integer           Number of CPUs to use per job ($opts{cpu})
  --job       integer           Number of simultaneous jobs ($opts{job})

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
  --meanbdg   path             ($opts{meanbdg})
  --bedtools  path             ($opts{bedtools})
  --intersect path             ($opts{intersect})
  --combrep   path             ($opts{combrep})
  --plotpeak  path             ($opts{plotpeak})
  --rscript   path             ($opts{rscript})
DOC


### Inputs
unless (@ARGV) {
	print $documentation;
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
	'chrskip=s'             => \$opts{chrskip},
	'blacklist=s'           => \$opts{blacklist},
	'dup!'                  => \$opts{dup},
	'dupfrac=f'             => \$opts{dupfrac},
	'maxdup=i'              => \$opts{maxdup},
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
	'bam2wig=s'             => \$opts{bam2wig},
	'bamdedup=s'            => \$opts{bamdedup},
	'macs=s'                => \$opts{macs},
	'manwig=s'              => \$opts{manwig},
	'wig2bw=s'              => \$opts{wig2bw},
	'bw2bdg=s'              => \$opts{bw2bdg},
	'bedtools=s'            => \$opts{bedtools},
	'getdata=s'             => \$opts{getdata},
	'printchr=s'            => \$opts{printchr},
	'meanbdg=s'             => \$opts{meanbdg},
	'intersect=s'           => \$opts{intersect},
) or die "unrecognized option(s)!\n";



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
run_dedup();
check_bams();
run_bam_conversion();
check_control();
run_bw_conversion();
run_bdgcmp();
run_call_peaks();
run_bdg_conversion();
run_peak_merge();
run_rescore();
run_plot_peaks();
finish();




############### Subroutines ########################################################

sub check_inputs {
	if (@ARGV) {
		die sprintf("There are unrecognized leftover items on the command line!\n %s\n", join("\n", @ARGV));
	}
	unless (@chips) {
		die "no ChIP file(s) defined!\n";
	}
	unless (@names) {
		die "no name(s) defined!\n";
	}
	unless (scalar(@chips) == scalar(@names)) {
		die "unequal ChIP and name quantities!\n";
	}
	if (scalar(@controls) > 1 and scalar(@controls) != scalar(@chips)) {
		die "Control and ChIP samples quantity doesn't match!\n";
	}
	elsif (scalar(@controls) == 0) {
		# no controls, turn off lambda
		$opts{slocal} = 0;
		$opts{llocal} = 0;
		$opts{lambda} = 0;
	}
	if (scalar(@chip_scales) and scalar(@chip_scales) != scalar(@chips)) {
		die "unequal ChIP samples and ChIP scale factors!\n";
	}
	if (scalar(@control_scales) and scalar(@control_scales) != scalar(@controls)) {
		die "unequal control samples and control scale factors!\n";
	}
	if (scalar(@chrnorms) and not $opts{chrapply}) {
		die "chromosome normalization factors given but no chromosome specified!\n";
	}
	if (not scalar(@chrnorms) and $opts{chrapply}) {
		die "chromosome name for normalization specified but no factors given!\n";
	}
	if (scalar(@chrnorms) and scalar(@chrnorms) != scalar(@chips)) {
		# apply to all the ChIPs
		if (scalar @chrnorms == 1) {
			print "WARNING: using the same chromosome normalization factor for each ChIP sample\n";
			my $n = shift @chrnorms;
			foreach (@names) {
				push @chrnorms, $n;
			}
		}
		else {
			die "ERROR: inequal number of chromosome normalization factors and ChIP samples!\n";
		}
	}
	if (scalar(@chrnorms) > 1 and scalar(@controls) == 1) {
		print "WARNING: using first chromosome normalization factor for universal control!\n";
	}
	if (not $opts{genome}) {
		my $s = $opts{species};
		$opts{genome} = $s eq 'human' ? 2700000000 : $s eq 'mouse' ? 1870000000 : 
			$s eq 'zebrafish' ? 1300000000 : $s eq 'fly' ? 120000000 : 
			$s eq 'celegens' ? 90000000 : $s eq 'yeast' ? 12100000 : 
			$s eq 'sheep' ? 2587000000 : 0;
		die "unknown species!\n" unless $opts{genome};
	}
	unless (-e $opts{dir} and -d _) {
		make_path($opts{dir}) or die sprintf("unable to make directory %s! $!\n", $opts{dir});
	}
	unless (-w $opts{dir}) {
		die sprintf("target directory %s is not writable!\n", $opts{dir});
	}
	unless (defined $opts{peaksize}) {
		$opts{peaksize} = 2 * $opts{fragsize};
	}
	unless (defined $opts{peakgap}) {
		$opts{peakgap} = $opts{fragsize};
	}
	unless (defined $opts{broadgap}) {
		$opts{broadgap} = 4 * $opts{fragsize};
	}
	if ($opts{plot} and scalar(@names) == 1) {
		# no sense plotting figures if we only have one ChIP set
		$opts{plot} = 0;
	}
	# add parameters to option hash for printing configuration
	$opts{chipscale} = join(", ", @chip_scales);
	$opts{controlscale} = join(", ", @control_scales);
	$opts{chrnorms} = join(", ", @chrnorms);
}

sub print_start {
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
		foreach my $f (split ',', $chips[$i]) {
			unless (-e $f) {
				printf(" can't find %s ChIP file %s!\n", $names[$i], $f);
				$error++;
			}
		}
		if ($controls[$i]) {
			foreach my $f (split ',', $controls[$i]) {
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
		my $universal_scale = shift @control_scales || undef;
		my $universal_name = 'Custom-Universal-' . $universal_control;
		foreach (@names) {
			push @controls, $universal_name;
		}
		push @jobs, ChIPjob->new($opts{out} . "_control", '', $universal_control, 
			undef, $universal_scale, $chrnorms[0] || undef);
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
		deduplication => 0,
		bam2wig       => 0,
		lambda        => 0,
		bw2bdg        => 0,
		bdgcmp        => 0,
		callpeak      => 0,
		bdg2bw        => 0,
		peakmerge     => 0,
		rescore       => 0
	);
	if (-e $progress_file) {
		my $fh = IO::File->new($progress_file, '<');
		while (my $line = $fh->getline) {
			chomp $line;
			$p{$line} = 1 if exists $p{$line};
		}
		$fh->close;
	}
	return %p;
}

sub run_dedup {
	return unless ($opts{dup});
	my @commands;
	my %name2done;
	print "\n\n======= De-duplicating bam files\n";
	if ($progress{deduplication}) {
		print "\nStep is completed\n";
		return;
	}
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_dedup_commands(\%name2done);
	}
	if (@commands) {
		execute_commands(\@commands);
	}
	update_progress_file('deduplication');
}

sub execute_commands {
	my $commands = shift;
	printf "Excecuting %d commands\n", scalar @$commands;
	if ($parallel) {
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
			print "=== Job: $command\n\n";
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
	
	# command
	my ($command_string, $command_out, $command_log) = @$command;
	my $command_app;
	if ($command_string =~ m/^([\w\_\.\/]+) /) {
		$command_app = $1;
	}
	
	# check
	if (length($command_out) and length($command_log)) {
		# both 
		if (-e $command_out and -e $command_log) {
			print "=== Job: $command_app previously finished, have $command_out and $command_log files\n" if $talk;
			return 1;
		}
		elsif (not -e $command_out and -e $command_log and 
			$command_app eq $opts{bamdedup}) 
		{
			# the deduplication command will not write out a bam file if the actual 
			# duplication rate is below the target rate
			# presume this is good!?
			print "=== Job: $command_app presumed finished, have $command_log file only\n" if $talk;
			return 2;
		}
		elsif (-e $command_out and not -e $command_log) {
			# we have a output file but not a log file
			print "=== Job: $command_app presumed finished, have $command_out file only, no log\n" if $talk;
			return 3;
		}
	}
	elsif (length($command_out)) {
		if (-e $command_out) {
			print "=== Job: $command_app previously finished, have $command_out\n" if $talk;
			return 4;
		}
	}
	elsif ($command_app eq 'rm') {
		# remove command doesn't leave an output (duh!) or log file
		# gotta check each one
		my $check = 0; 
		foreach my $item (split /\s+/, $command_string) {
			next if $item eq 'rm';
			$check++ if -e $item; # check true if file is present
		}
		if ($check == 0) {
			print "=== Job: $command_app previously finished, target files missing\n" if $talk;
			return 5;
		}
	}
	return 0; # command was not finished
}

sub update_progress_file {
	my $key = shift;
	$progress{$key} = 1;
	my $fh = IO::File->new($progress_file, '>>') or 
		die "can't write to progress file! $!\n";
	$fh->print("$key\n");
	$fh->close;
	return 1;
}

sub check_bams {
	print "\n\n======= Checking bam files\n";
	foreach my $Job (@Jobs) {
		$Job->find_dedup_bams;
	}
}

sub run_bam_conversion {
	my @commands;
	my %name2done;
	print "\n\n======= Converting bam files\n";
	if ($progress{bam2wig}) {
		print "\nStep is completed\n";
		return;
	}
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_bam2wig_commands(\%name2done);
	}
	if (@commands) {
		execute_commands(\@commands);
	}
	update_progress_file('bam2wig');
}

sub check_control {
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

sub generate_chr_file {
	my @bams = split(',', $chips[0]);
	my $example = shift @bams;
		# this will work regardless if example is bam or bigWig
	my $chromofile = File::Spec->catfile($opts{dir},"chrom_sizes.temp.txt");
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
	my $command = sprintf("%s --bed %s --out %s ", $opts{intersect}, $opts{bedtools}, 
		 $merge_file);
	foreach my $Job (@Jobs) {
		if ($Job->{peak}) {
			$command .= sprintf("%s ", $Job->{peak});
		}
	}
	my $log = $merge_file . '.merge.out.txt';
	$command .= sprintf("2>&1 > %s ", $log);
	push @commands, [$command, $merge_file . '.bed', $log]; 
		# this will have multiple outputs, but one is just a .bed file
	
	# broadPeaks
	if ($opts{broad}) {
		my $merge2_file = File::Spec->catfile($opts{dir}, $opts{out} . "_broad");
		my $command2 = sprintf("%s --bed %s --out %s ", $opts{intersect}, $opts{bedtools}, 
			 $merge2_file);
	
		foreach my $Job (@Jobs) {
			if ($Job->{peak}) {
				my $bf = $Job->{peak};
				$bf =~ s/narrow/gapped/;
				$command2 .= sprintf("%s ", $bf);
			}
		}
		my $log2 = $merge2_file . '.merge.out.txt';
		$command2 .= sprintf("2>&1 > %s ", $log2);
		push @commands, [$command2, $merge2_file . '.bed', $log2];
			# this will have multiple outputs, but one is just a .bed file
	}
	
	execute_commands(\@commands);
	update_progress_file('peakmerge');
}

sub run_rescore {
	print "\n\n======= Re-scoring all merged peaks\n";
	if ($progress{rescore}) {
		print "\nStep is completed\n";
		return;
	}
	die "no get_datasets.pl script in path!\n" unless $opts{getdata} =~ /\w+/;
	my %name2done;
	
	# prepare filenames
	my $input;
	if (scalar(@Jobs) > 1) {
		$input = File::Spec->catfile($opts{dir}, $opts{out} . '.bed');
		die "unable to find merged bed file '$input'!\n" unless -e $input;
	}
	else {
		$input = $Jobs[0]->{peak};
		die "unable to find peak file '$input'!\n" unless -e $input;
	}
	my $output1 = File::Spec->catfile($opts{dir}, $opts{out} . '_qvalue.txt');
	my $output2 = File::Spec->catfile($opts{dir}, $opts{out} . '_log2FE.txt');
	my $output3 = File::Spec->catfile($opts{dir}, $opts{out} . '_counts.txt');
	my $output4 = File::Spec->catfile($opts{dir}, $opts{out} . '_genome_counts.txt.gz');
	
	# start list of conditions
	my @conditions = ("Sample\tGroup\n");
	
	# generate four get_dataset commands
	# go ahead and make the fourth genome-wide command, even though we may not use it
	my $command1 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
		$opts{getdata}, $opts{cpu}, $input, $output1);
	my $command2 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
		$opts{getdata}, $opts{cpu}, $input, $output2);
	my $command3 = sprintf("%s --method sum --cpu %s --in %s --out %s --format 0 ",
		$opts{getdata}, $opts{cpu}, $input, $output3);
	my $command4 = sprintf("%s --method sum --cpu %s --feature genome --win %d --discard %s --out %s --format 0 ",
		$opts{getdata}, $opts{cpu}, $opts{window}, $opts{discard}, $output4);
	foreach my $Job (@Jobs) {
		if ($Job->{qvalue_bw}) {
			$command1 .= sprintf("--data %s ", $Job->{qvalue_bw});
		}
		if ($Job->{logfe_bw}) {
			$command2 .= sprintf("--data %s ", $Job->{logfe_bw});
		}
		foreach my $b ( @{ $Job->{chip_count_bw} } ) {
			$command3 .=  "--data $b ";
			$command4 .=  "--data $b ";
			# get the sample name just as it would appear in get_datasets output... a hack
			my (undef, undef, $name) = File::Spec->splitpath($b);
			$name =~ s/^([\w\d\-\_]+)\..+$/$1/i; # take everything up to first .
			push @conditions, sprintf("%s\t%s\n", $name, $Job->{name})
		}
		foreach my $b ( @{ $Job->{control_count_bw} } ) {
			next if exists $name2done{$b};
			$command3 .=  "--data $b ";
			$command4 .=  "--data $b ";
			my (undef, undef, $name) = File::Spec->splitpath($b);
			$name =~ s/^([\w\d\-\_]+)\..+$/$1/i; # take everything up to first .
			push @conditions, "$name\tInput\n";
			$name2done{$b} = 1; # remember it's done
		}
	}
	
	# add log outputs to commands
	my @commands;
	my $log = $output1;
	$log =~ s/txt$/log.txt/;
	$command1 .= " 2>&1 > $log";
	push @commands, [$command1, $output1, $log];
	$log = $output2;
	$log =~ s/txt$/log.txt/;
	$command2 .= " 2>&1 > $log";
	push @commands, [$command2, $output2, $log];
	$log = $output3;
	$log =~ s/txt$/log.txt/;
	$command3 .= " 2>&1 > $log";
	push @commands, [$command3, $output3, $log];
	if ($opts{window}) {
		# user has given an actual genome window size, so we'll run this command
		$log = $output4;
		$log =~ s/txt\.gz$/log.txt/;
		$command4 .= " 2>&1 > $log";
		push @commands, [$command4, $output4, $log];
	}
	
	# broad peak rescore
	my ($output5, $output6, $output7);
	if ($opts{broad}) {
		# generate broad file names
		my $input2;
		if (scalar(@Jobs) > 1) {
			$input2 = File::Spec->catfile($opts{dir}, $opts{out} . '_broad.bed');
			die "unable to find merged bed file '$input2'!\n" unless -e $input2;
		}
		else {
			$input2 = $Jobs[0]->{peak};
			$input2 =~ s/narrow/gapped/;
			die "unable to find merged bed file '$input2'!\n" unless -e $input2;
		}
		$output5 = File::Spec->catfile($opts{dir}, $opts{out} . '_broad_qvalue.txt');
		$output5 = File::Spec->catfile($opts{dir}, $opts{out} . '_broad_log2FE.txt');
		$output7 = File::Spec->catfile($opts{dir}, $opts{out} . '_broad_counts.txt');
		
		# generate three get_dataset commands
		my $command5 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
			$opts{getdata}, $opts{cpu}, $input2, $output5);
		my $command6 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
			$opts{getdata}, $opts{cpu}, $input2, $output6);
		my $command7 = sprintf("%s --method sum --cpu %s --in %s --out %s --format 0 ",
			$opts{getdata}, $opts{cpu}, $input2, $output7);
		%name2done = ();
		foreach my $Job (@Jobs) {
			if ($Job->{qvalue_bw}) {
				$command5 .= sprintf("--data %s ", $Job->{qvalue_bw});
			}
			if ($Job->{logfe_bw}) {
				$command6 .= sprintf("--data %s ", $Job->{logfe_bw});
			}
			foreach my $b ( @{ $Job->{chip_count_bw} } ) {
				$command7 .=  "--data $b ";
			}
			foreach my $b ( @{ $Job->{control_count_bw} } ) {
				next if exists $name2done{$b};
				$command7 .=  "--data $b ";
				$name2done{$b} = 1; # remember it's done
			}
		}
	
		# add log outputs to commands
		my @commands;
		my $log = $output5;
		$log =~ s/txt$/log.txt/;
		$command5 .= " 2>&1 > $log";
		push @commands, [$command5, $output5, $log];
		$log = $output6;
		$log =~ s/txt$/log.txt/;
		$command6 .= " 2>&1 > $log";
		push @commands, [$command6, $output6, $log];
		$log = $output7;
		$log =~ s/txt$/log.txt/;
		$command7 .= " 2>&1 > $log";
		push @commands, [$command7, $output7, $log];
	}
	
	# write conditions file
	my $output8 = File::Spec->catfile($opts{dir}, $opts{out} . '_samples.txt');
	my $fh = IO::File->new($output8, "w");
	foreach (@conditions) {
		$fh->print($_);
	}
	$fh->close;
	
	
	
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
		my $command9 = sprintf("%s --in %s --out %s --sample %s --method mean --format 0 ", 
			$opts{combrep}, $output3, $output3m, $output8);
		$log = $output3m;
		$log =~ s/txt/log.txt/;
		$command9 .= " 2>&1 > $log";
		push @commands2, [$command9, $output3m, $log];
		
		# genome counts
		if ($opts{window}) {
			my $output4m = File::Spec->catfile($opts{dir}, $opts{out} . 
				'_genome_meanCounts.txt.gz');
			my $command10 = sprintf("%s --in %s --out %s --sample %s --method mean --format 0 ", 
				$opts{combrep}, $output4, $output4m, $output8);
			$log = $output4m;
			$log =~ s/txt\.gz/log.txt/;
			$command10 .= " 2>&1 > $log";
			push @commands2, [$command10, $output4m, $log];
		}
		
		# broadPeak Counts
		if ($opts{broad}) {
			my $output6m = File::Spec->catfile($opts{dir}, $opts{out} . 
				'_broad_meanCounts.txt');
			my $command11 = sprintf("%s --in %s --out %s --sample %s --method mean --format 0 ", 
				$opts{combrep}, $output6, $output6m, $output8);
			$log = $output6m;
			$log =~ s/txt/log.txt/;
			$command11 .= " 2>&1 > $log";
			push @commands2, [$command11, $output6m, $log];
		}
		
		# execute
		execute_commands(\@commands2);
	}
	update_progress_file('rescore');
}


sub run_plot_peaks {
	return unless ($opts{plot});
	print "\n\n======= Plotting Peak figures\n";
	if ($opts{rscript} !~ /\w+/ or $opts{plotpeak} !~ /\w+/) {
		print "Rscript or plot_peak_figures.R script not defined!\n";
		return;
	}
	my $outbase = File::Spec->catfile($opts{dir}, $opts{out});
	my $command = sprintf("%s %s --input %s ", $opts{rscript}, $opts{plotpeak}, $outbase);
	my $log = $outbase . '_plot_figures.out.txt';
	$command .= " 2>&1 > $log";
	execute_commands( [ [$command, '', $log] ] );
}


sub finish {
	# combine output logs
	my @combined_output;
	foreach my $c (@finished_commands) {
		my $log = $c->[2]; # the log file
		next unless $log;
		push @combined_output, "=== Log file: $log\n";
		if (-e $log) {
			if (-z _) {
				# an empty file
				push @combined_output, "\n";
				unlink $log;
			}
			else {
				# push log contents to combined output
				my $fh = IO::File->new($log) or next;
				push @combined_output, <$fh>;
				push @combined_output, "\n";
				$fh->close;
			}
			unlink $log if -e $log;
		}
	}
	my $file = File::Spec->catfile($opts{dir}, $opts{out} . "_job_output_logs.txt");
	my $fh = IO::File->new($file, "w");
	foreach (@combined_output) {
		$fh->print($_);
	}
	$fh->close;
	
	# remove files no longer need
	unlink $chromofile;
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
	unlink $progress_file;
	
	# final print statements
	print "\nCombined all job output log files into '$file'\n";
	printf "Finished in %.1f minutes\n", (time -$start) / 60;
	print "\n======== Finished ChIPSeq multi-replicate pipeline ==========\n\n";
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
	    name => $name,
	    chip_bams => undef,
	    control_bams => undef,
	    chip_dedup_bams => [],
	    chip_use_bams => [],
	    chip_scale => undef,
	    control_dedup_bams => [],
	    control_use_bams => [],
	    control_scale => undef,
	    chip_bw => undef,
	    chip_bdg => undef,
	    chip_count_bw => [],
	    d_control_bdg => undef,
	    control_count_bw => [],
	    s_control_bdg => undef,
	    l_control_bdg => undef,
	    sl_control_bdg => undef,
	    sld_control_file => undef,
	    lambda_bdg => undef,
	    lambda_bw => undef,
	    qvalue_bdg => undef,
	    fe_bdg => undef,
	    logfe_bw => undef,
	    peak => undef,
	    qvalue_bw => undef,
	    chrnorm => $chrnorm,
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
		if ($chip_scale) {
			$self->{chip_scale} = [split(',', $chip_scale)];
			$self->crash("unequal scale factors and bam files!\n") if 
				scalar(@{$self->{chip_bams}}) != scalar(@bams);
		}
		# count bw files
		foreach my $bam (@bams) {
			my $bamname = $bam; # $bam is aliased, so must make a copy
			$bamname =~ s/\.bam$//i;
			my (undef, undef, $fname) = File::Spec->splitpath($bamname);
			push @{ $self->{chip_count_bw} }, sprintf("%s.count.bw",
				File::Spec->catfile($opts{dir}, $fname));
		}
	}
	else {
		# must be a control
		$self->{chip_bams} = [];
	}
	
	
	## check the control files
	if ($control =~ /^Custom\-Universal\-(.+)$/) {
		# this is just a ChIP Job using a common universal reference file
		my (undef, undef, $file) = File::Spec->splitpath($1);
		$file =~ s/(?:bw|bigwig)$/bdg/i;
		$self->{lambda_bdg} = File::Spec->catfile($opts{dir}, $file); 
	}
	elsif ($control =~ /\.(?:bw|bigwig)$/i) {
		# pre-generated reference bigWig file
		$self->crash("only one control bigWig file is allowed per experiment!\n") if
			$chip =~ /,/;
		$self->{lambda_bw} = $control;
		my (undef, undef, $file) = File::Spec->splitpath($control);
		$file =~ s/(?:bw|bigwig)$/bdg/i;
		$self->{lambda_bdg} = File::Spec->catfile($opts{dir}, $file); 
	}
	elsif ($control =~ /\.bam$/i) {
		# we have control bams to process
		my @bams = split(',', $control);
		$self->{control_bams} = \@bams;
		my $lambdabase = $name =~ /_control$/ ? $namepath : $namepath . "_control";
			# make sure we get _control in the name
			# universal controls will get it, but individual controls will not
		if ($opts{lambda}) {
			$self->{lambda_bdg} = "$lambdabase.lambda_control.bdg";
			$self->{lambda_bw} = "$lambdabase.lambda_control.bw";
			$self->{d_control_bdg} = "$lambdabase.dlocal.bdg";
			$self->{s_control_bdg} = "$lambdabase.slocal.bdg" if $opts{slocal};
			$self->{l_control_bdg} = "$lambdabase.llocal.bdg" if $opts{llocal};
			$self->{sl_control_bdg} = "$lambdabase.sllocal.bdg";
			$self->{sld_control_file} = "$lambdabase.sldlocal.txt";
		}
		else {
			$self->{lambda_bdg} = $lambdabase . ".bdg";
			$self->{lambda_bw} = $lambdabase . ".bw";
		}
		
		if ($control_scale) {
			$self->{control_scale} = [split(',', $control_scale)];
			$self->crash("unequal scale factors and bam files!\n") if 
				scalar(@bams) != scalar(@{$self->{control_scale}});
		}
		# count bw files
		foreach my $bam (@bams) {
			my $bamname = $bam; # $bam is aliased, so must make a copy
			$bamname =~ s/\.bam$//i;
			my (undef, undef, $fname) = File::Spec->splitpath($bamname);
			push @{ $self->{control_count_bw} }, sprintf("%s.count.bw",
				File::Spec->catfile($opts{dir}, $fname));
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
			my (undef,undef,$out) = File::Spec->splitpath($in);
			$out = File::Spec->catfile($opts{dir}, $out);
			$out =~ s/\.bam$/.dedup.bam/i;
			$self->{chip_dedup_bams}->[$i] = $out;
			my $command = sprintf "%s --in %s --out %s --cpu %s ", 
				$opts{bamdedup},
				$in,
				$out,
				$opts{cpu};
			if ($opts{dupfrac} > 0) {
				$command .= sprintf("--random --seed 1 --frac %s ", $opts{dupfrac});
			}
			if (defined $opts{maxdup}) {
				$command .= sprintf("--max %s ", $opts{maxdup});
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
			my (undef,undef,$out) = File::Spec->splitpath($in);
			$out = File::Spec->catfile($opts{dir}, $out);
			$out =~ s/\.bam$/.dedup.bam/i;
			$self->{control_dedup_bams}->[$i] = $out;
			my $command = sprintf "%s --in %s --out %s --cpu %s ", 
				$opts{bamdedup},
				$in,
				$out,
				$opts{cpu};
			if ($opts{dupfrac} > 0) {
				$command .= sprintf("--random --seed 1 --frac %s ", $opts{dupfrac});
			}
			if (defined $opts{maxdup}) {
				$command .= sprintf("--max %s ", $opts{maxdup});
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

sub generate_bam2wig_commands {
	my $self = shift;
	my $name2done = shift;
	croak("no bam2wig.pl script in path!\n") unless $opts{bam2wig} =~ /\w+/;
	my @commands;
	
	### ChIP bams
	if (scalar @{$self->{chip_use_bams}}) {
		# we have bam files to convert to bw
		
		# base fragment command
		my $frag_command = sprintf(
			"%s --in %s --out %s --qual %s --nosecondary --noduplicate --nosupplementary --bin %s --cpu %s --rpm --mean --bw --bwapp %s ", 
			$opts{bam2wig}, 
			join(',', @{$self->{chip_use_bams}}), 
			$self->{chip_bw},
			$opts{mapq},
			$opts{chipbin},
			$opts{cpu},
			$opts{wig2bw}
		);
		
		# base count command
		my $count_command = sprintf(
			"%s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --bw --bwapp %s ", 
			$opts{bam2wig}, 
			$opts{mapq},
			$opts{cpu},
			$opts{wig2bw}
		);
		# paired or single options
		if ($opts{paired}) {
			$frag_command .= sprintf("--pe --span --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
			$count_command .= sprintf("--pe --mid --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
		}
		else {
			$frag_command .= sprintf("--extend --extval %s ", $opts{fragsize});
			$frag_command .= sprintf("--shiftval %s ", $opts{shiftsize}) 
				if $opts{shiftsize};
			$count_command .= sprintf("--start --shiftval %0.0f ", 
				($opts{fragsize} / 2) + $opts{shiftsize} );
		}
		# additional filtering
		if ($opts{blacklist}) {
			$frag_command .= sprintf("--blacklist %s ", $opts{blacklist});
			$count_command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$frag_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
			$count_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		# scaling
		if ($self->{chip_scale}) {
			$frag_command .= sprintf("--scale %s ", join(',', @{$self->{chip_scale}} ) );
			# count command gets scaled below
		}
		# chromosome-specific scaling
		if ($self->{chrnorm}) {
			$frag_command .= sprintf("--chrnorm %s --chrapply %s ", $self->{chrnorm}, 
				$opts{chrapply});
			# count command gets scaled below
		}
		# finish fragment command
		my $log = $self->{chip_bw};
		$log =~ s/bw$/out.txt/;
		$frag_command .= " 2>&1 > $log";
		push @commands, [$frag_command, $self->{chip_bw}, $log];
		
		# finish count commands
		for my $i (0 .. $#{$self->{chip_use_bams}} ) {
			# add filenames
			my $command = $count_command . sprintf("--in %s --out %s ", 	
				$self->{chip_use_bams}->[$i], $self->{chip_count_bw}->[$i]);
			# add scaling as necessary
			unless ($opts{rawcounts}) {
				$command .= "--rpm --format 3 ";
				if ($self->{chrnorm}) {
					# chromosome-specific scaling
					$count_command .= sprintf("--chrnorm %s --chrapply %s ", 
						$self->{chrnorm}, $opts{chrapply});
				}
				if ($self->{chip_scale}) {
					$command .= sprintf("--scale %.4f ", 
						$opts{targetdep} * $self->{chip_scale}->[$i] );
				}
				else {
					# we always scale the count by the target depth
					$command .= sprintf("--scale %d ", $opts{targetdep})
				}
			}
			my $log = $self->{chip_count_bw}->[$i];
			$log =~ s/bw$/out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [$command, $self->{chip_count_bw}->[$i], $log];
		}
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
			"%s --in %s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --rpm --mean --bdg ", 
			$opts{bam2wig}, 
			$control_bam_string, 
			$opts{mapq},
			$opts{cpu}
		);
		# base count command
		my $count_command = sprintf(
			"%s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --bw --bwapp %s ", 
			$opts{bam2wig}, 
			$opts{mapq},
			$opts{cpu},
			$opts{wig2bw}
		);
		
		# general paired options, restrict size for all
		if ($opts{paired}) {
			$frag_command .= sprintf("--pe --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
			$count_command .= sprintf("--pe --mid --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
		}
		else {
			$count_command .= sprintf("--start --shiftval %0.0f ", 
				($opts{fragsize} / 2) + $opts{shiftsize});
		}
		
		# additional filters
		if ($opts{blacklist}) {
			$frag_command .= sprintf("--blacklist %s ", $opts{blacklist});
			$count_command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$frag_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
			$count_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		# chromosome-specific scaling
		if ($self->{chrnorm}) {
			$frag_command .= sprintf("--chrnorm %s --chrapply %s ", $self->{chrnorm}, 
				$opts{chrapply});
			# count command gets scaled below
		}
		
		# add count commands
		for my $i (0 .. $#{$self->{control_use_bams}} ) {
			# add filenames
			my $command = $count_command . sprintf("--in %s --out %s ", 	
				$self->{control_use_bams}->[$i], $self->{control_count_bw}->[$i]);
			# add scaling as necessary
			unless ($opts{rawcounts}) {
				$command .= "--rpm --format 3 ";
				if ($self->{chrnorm}) {
					# chromosome-specific scaling
					$count_command .= sprintf("--chrnorm %s --chrapply %s ", 
						$self->{chrnorm}, $opts{chrapply});
				}
				if ($self->{control_scale}) {
					$command .= sprintf("--scale %.4f ", 
						$opts{targetdep} * $self->{control_scale}->[$i] );
				}
				else {
					# we always scale the count by the target depth
					$command .= sprintf("--scale %d ", $opts{targetdep})
				}
			}
			my $log = $self->{control_count_bw}->[$i];
			$log =~ s/bw$/out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [$command, $self->{control_count_bw}->[$i], $log];
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
		$command1 .= sprintf("--bin %s --out %s ", $opts{chipbin}, $self->{d_control_bdg});
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
			$command2 .= sprintf("--cspan --extval %s --scale %s --bin %s --out %s ",
				$opts{slocal}, $scale, $opts{slocalbin}, 
				$self->{s_control_bdg});
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
			$command3 .= sprintf("--cspan --extval %s --scale %s --bin %s --out %s ",
				$opts{llocal}, $scale, $opts{llocalbin}, 
				$self->{l_control_bdg});
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
			"%s --in %s --out %s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --rpm --mean --bdg ", 
			$opts{bam2wig}, 
			$control_bam_string, 
			$self->{lambda_bdg},
			$opts{mapq},
			$opts{cpu}
		);
		
		# count command
		my $count_command = sprintf(
			"%s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --bw --bwapp %s ", 
			$opts{bam2wig}, 
			$opts{mapq},
			$opts{cpu},
			$opts{wig2bw}
		);
		
		# single or paired options
		if ($opts{paired}) {
			$frag_command .= sprintf("--span --pe --minsize %s --maxsize %s --bin %s ", 
				$opts{minsize}, $opts{maxsize}, $opts{chipbin});
			$count_command .= sprintf("--mid --pe --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
		}
		else {
			$frag_command .= sprintf("--extend --extval %s --shiftval %0.0f --bin %s ", 
				$opts{fragsize}, $opts{shiftsize}, $opts{chipbin});
			$count_command .= sprintf("--start --shiftval %0.0f ", 
				($opts{fragsize} / 2) + $opts{shiftsize});
		}
		
		# additional filters
		if ($opts{blacklist}) {
			$frag_command .= sprintf("--blacklist %s ", $opts{blacklist});
			$count_command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$frag_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
			$count_command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		# chromosome-specific scaling
		if ($self->{chrnorm}) {
			$frag_command .= sprintf("--chrnorm %s --chrapply %s ", $self->{chrnorm}, 
				$opts{chrapply});
			# count command gets scaled below
		}
		# scaling for the fragment command only, count below
		if ($self->{control_scale}) {
			$frag_command .= sprintf("--scale %s ", join(',', @{$self->{control_scale}} ) );
		}
		
		# add frag command
		my $log = $self->{lambda_bdg};
		$log =~ s/bdg$/out.txt/;
		$frag_command .= " 2>&1 > $log";
		push @commands, [$frag_command, $self->{lambda_bdg}, $log];
		
		# add count command
		for my $i (0 .. $#{$self->{control_use_bams}} ) {
			# add filenames
			my $command = $count_command . sprintf("--in %s --out %s ", 	
				$self->{control_use_bams}->[$i], $self->{control_count_bw}->[$i]);
			# add scaling as necessary
			unless ($opts{rawcounts}) {
				$command .= "--rpm --format 3 ";
				if ($self->{chrnorm}) {
					# chromosome-specific scaling
					$count_command .= sprintf("--chrnorm %s --chrapply %s ", 
						$self->{chrnorm}, $opts{chrapply});
				}
				if ($self->{control_scale}) {
					$command .= sprintf("--scale %.4f ", 
						$opts{targetdep} * $self->{control_scale}->[$i] );
				}
				else {
					# we always scale the count by the target depth
					$command .= sprintf("--scale %d ", $opts{targetdep})
				}
			}
			my $log = $self->{control_count_bw}->[$i];
			$log =~ s/bw$/out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [$command, $self->{control_count_bw}->[$i], $log];
		}
		
		# record that we've done this bam
		# we store the lambda bedgraph name because it might be reused again for another job
		$name2done->{$control_bam_string} = $self->{lambda_bdg};
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
		# this has to be done after the ChIP fragment bigWig conversion, since it 
		# enter a race condition if done with the other bam2wig jobs
		my $log = $self->{lambda_bdg};
		$log =~ s/bdg/out.txt/;
		my $command = sprintf("%s %s %s 2>&1 > $log", $opts{meanbdg}, $self->{chip_bw}, 
			$self->{lambda_bdg});
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
	$self->crash("no d control bedGraph file '$dfile'!\n") unless (-e $dfile);
	$self->crash("no small control bedGraph file '$sfile'!\n") if ($sfile and not -e $sfile);
	$self->crash("no large control bedGraph file '$lfile'!\n") if ($lfile and not -e $lfile);
	
	# generate background lambda bedgraph
	# go ahead and do this immediately because it's so simple
	my $background = sprintf("%.4f", ( 1_000_000 * $opts{fragsize} ) / $opts{genome} );
	my $background_bdg = $self->{lambda_bdg};
	$background_bdg =~ s/lambda_control/background/;
	my $infh = IO::File->new($chromofile); # use the chromosome file as source
	my $outfh = IO::File->new($background_bdg, "w");
	while (my $line = $infh->getline) {
		chomp $line;
		my ($chr, $end) = split /\s/, $line;
		$outfh->printf("%s\t0\t%s\t%s\n", $chr, $end, $background);
	}
	$infh->close;
	$outfh->close;
	
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
		$command .= sprintf("&& rm %s %s %s %s %s ", $dfile, $sfile, $lfile, 
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
		$command .= sprintf("&& rm %s %s %s %s ", $dfile, $sfile, 
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
		$command .= sprintf("&& rm %s %s %s %s ", $dfile, $lfile, 
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
	if ($self->{chip_bw} and -e $self->{chip_bw}) {
		my $log = $self->{chip_bdg};
		$log =~ s/bdg$/out.txt/;
		my $command = sprintf("%s %s %s 2>> $log", $opts{bw2bdg}, $self->{chip_bw}, 
			$self->{chip_bdg});
		push @commands, [$command, $self->{chip_bdg}, $log]
	}
	if ($self->{lambda_bw} and -e $self->{lambda_bw} and 
		not exists $name2done->{$self->{lambda_bdg}}
	) {
		my $log = $self->{lambda_bdg};
		$log =~ s/bdg$/out.txt/;
		my $command = sprintf("%s %s %s 2>> $log", $opts{bw2bdg}, $self->{lambda_bw}, 
			$self->{lambda_bdg});
		push @commands, [$command, $self->{lambda_bdg}, $log];
		$name2done->{$self->{lambda_bdg}} = 1; # mark as done
	}
	return @commands;
}

sub generate_enrichment_commands {
	my $self = shift;
	croak("no macs2 application in path!\n") unless $opts{macs} =~ /\w+/;
	my $chip = $self->{chip_bdg} || undef;
	my $lambda = $self->{lambda_bdg} || undef;
	return unless ($chip and $lambda);
	$self->crash("no ChIP fragment file $chip!\n") unless -e $chip;
	$self->crash("no control lambda fragment file $lambda!\n") unless -e $lambda;
	
	my $command = sprintf("%s bdgcmp -t %s -c %s -S %s -m qpois FE -o %s %s ", 
		$opts{macs}, $chip, $lambda, $opts{targetdep}, $self->{qvalue_bdg}, 
		$self->{fe_bdg});
	if (not $opts{lambda}) {
		$command .= "-p 1 "; # add a pseudo count of 1 when doing explicit comparisons
	}
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
	$self->crash("no qvalue bedGraph file $qtrack!\n") unless -e $qtrack;
	my $command = sprintf("%s bdgpeakcall -i %s -c %s -l %s -g %s --no-trackline -o %s ",
		$opts{macs}, $qtrack, $opts{cutoff}, $opts{peaksize}, $opts{peakgap}, 
		$self->{peak}
	);
	my $log = $self->{peak};
	$log =~ s/narrowPeak$/peakcall.out.txt/;
	$command .= " 2> $log";
	if ($opts{broad}) {
		# sneak this in as an extra command
		my $bpeak = $self->{peak};
		$bpeak =~ s/narrow/gapped/;
		my $command2 = sprintf("%s bdgbroadcall -i %s -c %s -C %s -l %s -g %s -G %s -o %s ",
			$opts{macs}, $qtrack, $opts{cutoff}, $opts{broadcut}, $opts{peaksize}, 
			$opts{peakgap}, $opts{broadgap}, $bpeak
		);
		my $log2 = $self->{peak};
		$log2 =~ s/narrowPeak$/broadcall.out.txt/;
		$command2 .= " 2> $log2";
		return ( [$command, $self->{peak}, $log], [$command2, $bpeak, $log2] );
	}
	return [$command, $self->{peak}, $log];
}

sub generate_bdg2bw_commands {
	my $self = shift;
	my $chromofile = shift;
	my $name2done = shift;
	
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
			croak("no wigToBigWig application in path!\n") unless $opts{wig2bw} =~ /\w+/;
			my $log = $self->{lambda_bw};
			$log =~ s/bw$/out.txt/;
			my $command = sprintf("%s %s %s %s 2>&1 > $log && rm %s", 
				$opts{wig2bw},
				$self->{lambda_bdg},
				$chromofile,
				$self->{lambda_bw},
				$self->{lambda_bdg},
			);
			push @commands, [$command, $self->{lambda_bw}, $log]
		}
		$name2done->{$self->{lambda_bdg}} = 1; # remember it's done
	}
	if ($self->{qvalue_bdg} and $self->{qvalue_bw}) {
		croak("no wigToBigWig application in path!\n") unless $opts{wig2bw} =~ /\w+/;
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
		croak("no wigToBigWig application in path!\n") unless $opts{wig2bw} =~ /\w+/;
		croak("no manipulate_wig.pl script in path!\n") unless $opts{manwig} =~ /\w+/;
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





