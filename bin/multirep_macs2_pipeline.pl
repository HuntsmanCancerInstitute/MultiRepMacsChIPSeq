#!/usr/bin/perl

use strict;
use IO::File;
use File::Spec;
use File::Which;
use Getopt::Long;

my $VERSION = 6.2;

my $parallel;
eval {
	require Parallel::ForkManager;
	$parallel = 1;
};
my $big_helper;
eval {
	require Bio::ToolBox::big_helper;
	Bio::ToolBox::big_helper->import('generate_chromosome_file');
	$big_helper = 1;
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
	dedup       => 1,
	maxdup      => undef,
	dupfrac     => 0.1,
	savebam     => 0,
	fragsize    => 300,
	shiftsize   => 0,
	slocal      => 1000,
	llocal      => 10000,
	qvalue      => 2,
	peaksize    => 300,
	peakgap     => 200,
	targetdep   => 25,
	use_lambda  => 1,
	chrskip     => "chrM|MT|lambda|Adapter|PhiX",
	blacklist   => undef,
	cpu         => 4,
	chipbin     => 10,
	slocalbin   => 50,
	llocalbin   => 100,
	chrapply    => undef,
	bam2wig     => sprintf("%s", which 'bam2wig.pl'),
	bamdedup    => sprintf("%s", which 'bam_partial_dedup.pl'),
	macs        => sprintf("%s", which 'macs2'),
	manwig      => sprintf("%s", which 'manipulate_wig.pl'),
	wig2bw      => sprintf("%s", which 'wigToBigWig'),
	bw2bdg      => sprintf("%s", which 'bigWigToBedGraph'),
	bedtools    => sprintf("%s", which 'bedtools'),
	getdata     => sprintf("%s", which 'get_datasets.pl'),
	printchr    => sprintf("%s", which 'print_chromosome_lengths.pl'),
	meanbdg     => sprintf("%s", which 'generate_mean_bedGraph.pl'),
	intersect   => sprintf("%s", which 'intersect_peaks.pl'),
) or die " unrecognized parameter!\n";
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
be provided as comma-delimited lists. If no control is available (for example, 
ATACSeq often has no genomic input), then a global mean coverage will be calculated 
from the ChIP samples and used as the control.

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
  --pe                          Bam files are paired-end
  --min       integer           Minimum paired-end size allowed ($opts{minsize} bp)
  --max       integer           Maximum paired-end size allowed ($opts{maxsize} bp)
 
 Bam filtering options
  --chrskip   "text"            Chromosome skip regex ($opts{chrskip})
  --blacklist file              Bed file of repeats or hotspots to avoid
  
 Duplication filtering
  --nodup                       Skip deduplication
  --dupfrac   fraction          Minimum allowed fraction of duplicates ($opts{dupfrac})
  --maxdup    integer           Maximum allowed duplication depth ($opts{maxdup})
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
  --cutoff    number            Threshold q-value for calling peaks ($opts{qvalue}) 
                                 Higher numbers are more significant, -1*log10(q)
  --tdep      integer           Average sequence depth of bam files in millions ($opts{targetdep})
  --peaksize  integer           Minimum peak size to call ($opts{peaksize} bp)
  --peakgap   integer           Maximum gap between peaks before merging ($opts{peakgap} bp)
  --nolambda                    Skip lambda control, compare ChIP directly with control
  
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
  --getdata   path             ($opts{getdata})
  --meanbdg   path             ($opts{meanbdg})
  --bedtools  path             ($opts{bedtools})
  --intersect path             ($opts{intersect})
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
	'pe!'                   => \$opts{paired},
	'min=i'                 => \$opts{minsize},
	'max=i'                 => \$opts{maxsize},
	'chrskip=s'             => \$opts{chrskip},
	'blacklist=s'           => \$opts{blacklist},
	'dup!'                  => \$opts{dedup},
	'dupfrac=f'             => \$opts{dupfrac},
	'maxdup=i'              => \$opts{maxdup},
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
	'cutoff=f'              => \$opts{qvalue},
	'tdep=f'                => \$opts{targetdep},
	'peaksize=i'            => \$opts{peaksize},
	'peakgap=i'             => \$opts{peakgap},
	'lambda!'               => \$opts{use_lambda},
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
my @finished_commands;
check_inputs();
print_start();
my $chromofile = generate_chr_file();
my @Jobs = generate_job_file_structure();
run_dedup() if ($opts{dedup});
check_bams();
run_bam_conversion();
check_control();
run_bw_conversion();
run_bdgcmp();
run_call_peaks();
run_bdg_conversion();
run_peak_merge();
run_rescore();
finish();




############### Subroutines ########################################################

sub check_inputs {
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
		mkdir $opts{dir} or die sprintf("unable to make directory %s! $!\n", $opts{dir});
	}
	unless (-w $opts{dir}) {
		die sprintf("target directory %s is not writable!\n", $opts{dir});
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

sub generate_job_file_structure {
	# we pass this off to ChIPjob package with 6 options: 
	# name, chip files, control files, chip scale factor, control scale factor, 
	# chromosome normalization factor
	my @jobs;
	# check for a unversal control
	if (scalar(@controls) == 1 and scalar(@chips) > 1) {
		print "Using only one control for multiple ChIP experiments\n";
		my $universal_control = shift @controls;
		my $universal_scale = shift @control_scales || undef;
		foreach (@names) {
			push @controls, 'universal';
		}
		push @jobs, ChIPjob->new($opts{out} . "_control", '', $universal_control, 
			undef, $universal_scale, $chrnorms[0] || undef);
	}
	for my $i (0 .. $#names) {
		push @jobs, ChIPjob->new($names[$i], $chips[$i], $controls[$i] || undef, 
			$chip_scales[$i] || undef, $control_scales[$i] || undef, 
			$chrnorms[$i] || undef );
	}
	return @jobs;
}

sub run_dedup {
	my @commands;
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_dedup_commands;
	}
	if (@commands) {
		print "\n\n======= De-duplicating bam files\n";
		execute_commands(\@commands);
	}
}

sub execute_commands {
	my $commands = shift;
	printf "Excecuting %d commands\n", scalar @$commands;
	if ($parallel) {
		my $pm = Parallel::ForkManager->new($opts{job});
		foreach my $command (@$commands) {
			next if check_command_finished($command);
			printf "=== Job: %s\n", $command->[0];
			$pm->start and next;
			# in child
			system($command->[0]);
			$pm->finish;
		}
		$pm->wait_all_children;
	}
	else {
		foreach my $command (@$commands) {
			next if check_command_finished($command);
			print "=== Job: $command\n\n";
			system($command->[0]);
		}
	}
	push @finished_commands, @$commands;
}

sub check_command_finished {
	my $command = shift;
	my ($command_string, $command_out, $command_log) = @$command;
	my $command_app;
	if ($command_string =~ m/^([\w\_\.\/]+) /) {
		$command_app = $1;
	}
	if (length($command_out) and length($command_log)) {
		# both 
		if (-e $command_out and -e $command_log) {
			print "=== Job: $command_app previously finished, have $command_out and $command_log files\n";
			return 1;
		}
		elsif (not -e $command_out and -e $command_log and 
			$command_app eq $opts{bamdedup}) 
		{
			# the deduplication command will not write out a bam file if the actual 
			# duplication rate is below the target rate
			# presume this is good!?
			print "=== Job: $command_app presumed finished, have $command_log file only\n";
			return 1;
		}
	}
	elsif (length($command_out)) {
		if (-e $command_out) {
			print "=== Job: $command_app previously finished, have $command_out\n";
			return 1;
		}
	}
	elsif ($command_app eq 'rm') {
		return; 
	}
	return;
}

sub check_bams {
	print "\n\n======= Checking bam files\n";
	foreach my $Job (@Jobs) {
		$Job->find_dedup_bams;
	}
}

sub run_bam_conversion {
	my @commands;
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_bam2wig_commands;
	}
	if (@commands) {
		print "\n\n======= Converting bam files\n";
		execute_commands(\@commands);
	}
}

sub check_control {
	return unless ($opts{use_lambda});
	my @commands;
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_lambda_control_commands;
	}
	if (@commands) {
		print "\n\n======= Generate control lambda files\n";
		execute_commands(\@commands);
	}
}

sub run_bw_conversion {
	my @commands;
	foreach my $Job (@Jobs) {
		push @commands, $Job->convert_bw_to_bdg;
	}
	if (@commands) {
		print "\n\n======= Converting Fragment bigWig files to bedGraph\n";
		execute_commands(\@commands);
	}
}

sub run_bdgcmp {
	my @commands;
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_enrichment_commands;
	}
	print "\n\n======= Generate enrichment files\n";
	execute_commands(\@commands);
}

sub run_call_peaks {
	my @commands;
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_peakcall_commands;
	}
	print "\n\n======= Call peaks\n";
	execute_commands(\@commands);
}

sub run_bdg_conversion {
	my @commands;
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_bdg2bw_commands($chromofile);
	}
	print "\n\n======= Converting bedGraph files\n";
	execute_commands(\@commands);
}

sub generate_chr_file {
	my $chromofile;
	my @bams = split(',', $chips[0]);
	my $example = shift @bams;
	if ($big_helper) {
		$chromofile = generate_chromosome_file($example, $opts{chrskip});
	}
	elsif ($opts{printchr}) {
		$chromofile = File::Spec->catfile($opts{dir},"chrom_sizes.temp.txt");
		system(sprintf("%s $example > $chromofile", $opts{printchr}));
		die "no chromosome file $chromofile!\n" unless -e $chromofile;
	}
	else {
		die "unable to generate chromosome file !\n";
	}
	return $chromofile;
}

sub run_peak_merge {
	return if scalar(@Jobs) == 1; # no sense merging one job!
	print "\n\n======= Merging called narrowPeak files\n";
	die "no bedtools application in path!\n" unless $opts{bedtools} =~ /\w+/;
	die "no intersect_peaks.pl application in path!\n" unless $opts{intersect} =~ /\w+/;
	my $merge_file = File::Spec->catfile($opts{dir}, $opts{out});
	my $command = sprintf("%s --bed %s --out %s ", $opts{intersect}, $opts{bedtools}, 
		 $merge_file);
	
	foreach my $Job (@Jobs) {
		if ($Job->{peak}) {
			$command .= sprintf("%s ", $Job->{peak});
		}
	}
	execute_commands([[$command, $merge_file, '']]);
}

sub run_rescore {
	print "\n\n======= Re-scoring all merged peaks\n";
	die "no get_datasets.pl script in path!\n" unless $opts{getdata} =~ /\w+/;
	
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
	
	# start list of conditions
	my @conditions = ("Sample\tCondition\n");
	
	# generate three get_dataset commands
	my $command1 = sprintf("%s --method mean --cpu %s --in %s --out %s ",
		$opts{getdata}, $opts{cpu}, $input, $output1);
	my $command2 = sprintf("%s --method mean --cpu %s --in %s --out %s ",
		$opts{getdata}, $opts{cpu}, $input, $output2);
	my $command3 = sprintf("%s --method sum --cpu %s --in %s --out %s ",
		$opts{getdata}, $opts{cpu}, $input, $output3);
	foreach my $Job (@Jobs) {
		if ($Job->{qvalue_bw}) {
			$command1 .= sprintf("--data %s ", $Job->{qvalue_bw});
		}
		if ($Job->{logfe_bw}) {
			$command2 .= sprintf("--data %s ", $Job->{logfe_bw});
		}
		foreach my $b ( @{ $Job->{chip_count_bw} } ) {
			$command3 .=  "--data $b ";
			# get the sample name just as it would appear in get_datasets output... a hack
			my (undef, undef, $name) = File::Spec->splitpath($b);
			$name =~ s/^([\w\d\-\_]+)\..+$/$1/i; # take everything up to first .
			push @conditions, sprintf("%s\t%s\n", $name, $Job->{name})
		}
		foreach my $b ( @{ $Job->{control_count_bw} } ) {
			$command3 .=  "--data $b ";
			my (undef, undef, $name) = File::Spec->splitpath($b);
			$name =~ s/^([\w\d\-\_]+)\..+$/$1/i; # take everything up to first .
			push @conditions, "$name\tInput\n";
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
	
	# write conditions file
	my $output4 = File::Spec->catfile($opts{dir}, $opts{out} . '_conditions.txt');
	my $fh = IO::File->new($output4, "w");
	foreach (@conditions) {
		$fh->print($_);
	}
	$fh->close;
	
	execute_commands(\@commands);
}

sub finish {
	# combine output logs
	my @combined_output;
	foreach my $c (@finished_commands) {
		my $log = $c->[2]; # the log file
		push @combined_output, "=== Log file: $log\n";
		if (-e $log) {
			if (-z _) {
				# an empty file
				push @combined_output, "\n";
				unlink $log
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
	
	# final print statements
	print "\nCombined all job output log files into '$file'\n";
	printf "Finished in %.1f minutes\n", (time -$start) / 60;
	print "\n======== Finished ChIPSeq multi-replicate pipeline ==========\n\n";
}



############### ChIPJob Package ########################################################

package ChIPjob;

sub new {
	# pass name, comma-list of ChIP files, and comma-list of control files
	my ($class, $name, $chip, $control, $chip_scale, $control_scale, $chrnorm) = @_;
	my $namepath = File::Spec->catfile($opts{dir}, $name);
	my $self = {
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
	    sld_control_bdg => undef,
	    lambda_bdg => undef,
	    lambda_bw => undef,
	    qvalue_bdg => undef,
	    fe_bdg => undef,
	    logfe_bw => undef,
	    peak => undef,
	    qvalue_bw => undef,
	    chrnorm => $chrnorm,
	};
	
	
	## check the ChIP files
	if ($chip =~ /\.(?:bw|bigwig)$/i) {
		die "only one ChIP bigWig file is allowed per experiment!\n" if $chip =~ /,/;
		$self->{chip_bw} = $chip;
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
			die "unequal scale factors and bam files!\n" if 
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
	if ($control eq 'universal') {
		# this is just a ChIP Job using a common universal control lambda
		$self->{lambda_bdg} = File::Spec->catfile($opts{dir}, 
			$opts{out} . "_control");
		$self->{lambda_bdg} .= $opts{use_lambda} ? '.lambda_control.bdg' : '.bdg';
	}
	elsif ($control =~ /\.(?:bw|bigwig)$/i) {
		die "only one control biwWig file is allowed per experiment!\n" if $chip =~ /,/;
		$self->{lambda_bw} = $control;
		$self->{lambda_bdg} = $control;
		$self->{lambda_bdg} =~ s/(?:bw|bigwig)$/bdg/i;
	}
	elsif ($control =~ /\.bam$/i) {
		# we have control bams to process
		my @bams = split(',', $control);
		$self->{control_bams} = \@bams;
		if ($opts{use_lambda}) {
			$self->{lambda_bdg} = "$namepath.lambda_control.bdg";
			$self->{lambda_bw} = "$namepath.lambda_control.bw";
			$self->{d_control_bdg} = "$namepath.dlocal.bdg";
			$self->{s_control_bdg} = "$namepath.slocal.bdg" if $opts{slocal};
			$self->{l_control_bdg} = "$namepath.llocal.bdg" if $opts{llocal};
			$self->{sl_control_bdg} = "$namepath.sllocal.bdg";
			$self->{sld_control_bdg} = "$namepath.sldlocal.bdg";
		}
		else {
			$self->{lambda_bdg} = "$namepath.bdg";
			$self->{lambda_bw} = "$namepath.bw";
		}
		if ($control_scale) {
			$self->{control_scale} = [split(',', $control_scale)];
			die "unequal scale factors and bam files!\n" if 
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
		$self->{lambda_bdg} = "$namepath.expected_mean.bdg";
		$self->{lambda_bw} = "$namepath.expected_mean.bw";
	}
	
	return bless $self, $class;
}

sub generate_dedup_commands {
	my $self = shift;
	die "no bam_partial_dedup.pl script in path!\n" unless $opts{bamdedup} =~ /\w+/;
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
			if ($opts{paired}) {
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
			if ($opts{paired}) {
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
	die "no bam2wig.pl script in path!\n" unless $opts{bam2wig} =~ /\w+/;
	my @commands;
	
	# ChIP bams
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
			"%s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --rpm --format 0 --bw --bwapp %s ", 
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
			$count_command .= sprintf("--chrnorm %s --chrapply %s ", $self->{chrnorm}, 
				$opts{chrapply});
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
			# add scaling
			if ($self->{chip_scale}) {
				$command .= sprintf("--scale %.4f ", 
					$opts{targetdep} * $self->{chip_scale}->[$i] );
			}
			else {
				# we always scale the count by the target depth
				$command .= sprintf("--scale %d ", $opts{targetdep})
			}
			my $log = $self->{chip_count_bw}->[$i];
			$log =~ s/bw$/out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [$command, $self->{chip_count_bw}->[$i], $log];
		}
	}
	
	# Control bams
	if (scalar @{$self->{control_use_bams}} and $opts{use_lambda}) {
		# process control bams into a chromatin bias lambda-control track
		
		# base fragment command
		my $frag_command = sprintf(
			"%s --in %s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --rpm --mean --bdg ", 
			$opts{bam2wig}, 
			join(',', @{$self->{control_use_bams}}), 
			$opts{mapq},
			$opts{cpu}
		);
		# count command
		my $count_command = sprintf(
			"%s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --rpm --mean --format 0 --bw --bwapp %s ", 
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
			$count_command .= sprintf("--chrnorm %s --chrapply %s ", $self->{chrnorm}, 
				$opts{chrapply});
		}
		
		# add count commands
		for my $i (0 .. $#{$self->{control_use_bams}} ) {
			# add filenames
			my $command = $count_command . sprintf("--in %s --out %s ", 	
				$self->{control_use_bams}->[$i], $self->{control_count_bw}->[$i]);
			# add scaling
			if ($self->{control_scale}) {
				$command .= sprintf("--scale %.4f ", 
					$opts{targetdep} * $self->{control_scale}->[$i] );
			}
			else {
				# we always scale the count by the target depth
				$command .= sprintf("--scale %d ", $opts{targetdep})
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
			$command1 .= sprintf("--extend --extval %s ", $opts{fragsize});
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
	}
	elsif (scalar @{$self->{control_use_bams}} and not $opts{use_lambda}) {
		# skipping chromatin-bias lambda control track, use control track as is
		
		# fragment command
		my $frag_command = sprintf(
			"%s --in %s --out %s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --rpm --mean --bdg ", 
			$opts{bam2wig}, 
			join(',', @{$self->{control_use_bams}}), 
			$self->{lambda_bdg},
			$opts{mapq},
			$opts{cpu}
		);
		
		# count command
		my $count_command = sprintf(
			"%s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --rpm --mean --format 0 --bw --bwapp %s ", 
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
			$frag_command .= sprintf("--extend --extval %s --bin %s ", $opts{fragsize}, 
				$opts{chipbin});
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
			$count_command .= sprintf("--chrnorm %s --chrapply %s ", $self->{chrnorm}, 
				$opts{chrapply});
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
			# add scaling
			if ($self->{control_scale}) {
				$command .= sprintf("--scale %.4f ", 
					$opts{targetdep} * $self->{control_scale}->[$i] );
			}
			else {
				# we always scale the count by the target depth
				$command .= sprintf("--scale %d ", $opts{targetdep})
			}
			my $log = $self->{control_count_bw}->[$i];
			$log =~ s/bw$/out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [$command, $self->{control_count_bw}->[$i], $log];
		}
	}
	
	# finished
	return @commands;
}

sub generate_lambda_control_commands {
	my $self = shift;
	die "no macs2 application in path!\n" unless $opts{macs} =~ /\w+/;
	
	my $dfile = $self->{d_control_bdg};
	my $sfile = $self->{s_control_bdg};
	my $lfile = $self->{l_control_bdg};
	return unless ($dfile); # controls with lambda will always have  d_control_bdg
		# ChIP jobs and non-lambda controls will not
	die "no d control bedGraph file '$dfile'!\n" unless (-e $dfile);
	die "no small control bedGraph file '$sfile'!\n" if ($sfile and not -e $sfile);
	die "no large control bedGraph file '$lfile'!\n" if ($lfile and not -e $lfile);
	
	my $log = $self->{lambda_bdg};
	$log =~ s/bdg$/out.txt/;
	my $command;
	
	if ($sfile and $lfile) {
		# first step
		$command = sprintf("%s bdgcmp -m max -t $sfile -c $lfile -o %s ", 
			$opts{macs}, $self->{sl_control_bdg});
		$command .= " 2> $log ";
	
		# second step
		$command .= sprintf("&& %s bdgcmp -m max -t %s -c $dfile -o %s ", 
			$opts{macs}, $self->{sl_control_bdg}, $self->{sld_control_bdg});
		$command .= " 2>> $log ";
	
		# third step
		my $background = ( 1_000_000 * $opts{fragsize} ) / $opts{genome};
		$command .= sprintf("&& %s bdgopt -m max -p $background -i %s -o %s ", 
			$opts{macs}, $self->{sld_control_bdg}, $self->{lambda_bdg});
		$command .= " 2>> $log ";
	
		# clean up
		$command .= sprintf("&& rm %s %s %s %s %s ", $dfile, $sfile, $lfile, 
			$self->{sl_control_bdg}, $self->{sld_control_bdg});
	}
	elsif ($sfile and not $lfile) {
		# first step
		$command = sprintf("%s bdgcmp -m max -t $sfile -c $dfile -o %s ", 
			$opts{macs}, $self->{sld_control_bdg});
		$command .= " 2> $log ";
	
		# second step
		my $background = ( 1_000_000 * $opts{fragsize} ) / $opts{genome};
		$command .= sprintf("&& %s bdgopt -m max -p $background -i %s -o %s ", 
			$opts{macs}, $self->{sld_control_bdg}, $self->{lambda_bdg});
		$command .= " 2>> $log ";
	
		# clean up
		$command .= sprintf("&& rm %s %s %s ", $dfile, $sfile, $self->{sld_control_bdg});
	}
	elsif (not $sfile and $lfile) {
		# first step
		$command = sprintf("%s bdgcmp -m max -t $lfile -c $dfile -o %s ", 
			$opts{macs}, $self->{sld_control_bdg});
		$command .= " 2> $log ";
	
		# second step
		my $background = ( 1_000_000 * $opts{fragsize} ) / $opts{genome};
		$command .= sprintf("&& %s bdgopt -m max -p $background -i %s -o %s ", 
			$opts{macs}, $self->{sld_control_bdg}, $self->{lambda_bdg});
		$command .= " 2>> $log ";
	
		# clean up
		$command .= sprintf("&& rm %s %s %s ", $dfile, $lfile, $self->{sld_control_bdg});
	}
	else {
		die "programming error! how did we get here with no sfile and no lfile????";
	}
	
	return [$command, $self->{lambda_bdg}, $log];
}

sub convert_bw_to_bdg {
	my $self = shift;
	die "no bigWigToBedGraph application in path!\n" unless $opts{bw2bdg} =~ /\w+/;
	my @commands;
	if ($self->{chip_bw} and -e $self->{chip_bw}) {
		my $log = $self->{chip_bdg};
		$log =~ s/bdg$/out.txt/;
		my $command = sprintf("%s %s %s 2>> $log", $opts{bw2bdg}, $self->{chip_bw}, 
			$self->{chip_bdg});
		push @commands, [$command, $self->{chip_bdg}, $log]
	}
	if ($self->{lambda_bw} and -e $self->{lambda_bw}) {
		my $log = $self->{lambda_bdg};
		$log =~ s/bdg$/out.txt/;
		my $command = sprintf("%s %s %s 2>> $log", $opts{bw2bdg}, $self->{lambda_bw}, 
			$self->{lambda_bdg});
		push @commands, [$command, $self->{lambda_bdg}, $log];
	}
	return @commands;
}

sub generate_enrichment_commands {
	my $self = shift;
	die "no macs2 application in path!\n" unless $opts{macs} =~ /\w+/;
	my $chip = $self->{chip_bdg} || undef;
	my $lambda = $self->{lambda_bdg} || undef;
	return unless ($chip and $lambda);
	die "no ChIP fragment file $chip!\n" unless -e $chip;
	die "no control lambda fragment file $lambda!\n" unless -e $lambda;
	
	my $command = sprintf("%s bdgcmp -t %s -c %s -S %s -m qpois FE -o %s %s ", 
		$opts{macs}, $chip, $lambda, $opts{targetdep}, $self->{qvalue_bdg}, 
		$self->{fe_bdg});
	if (not $opts{use_lambda}) {
		$command .= "-p 1 "; # add a pseudo count of 1 when doing explicit comparisons
	}
	my $log = $self->{qvalue_bdg};
	$log =~ s/bdg$/out.txt/;
	$command .= " 2> $log ";
	return [$command, $self->{qvalue_bdg}, $log];
}

sub generate_peakcall_commands {
	my $self = shift;
	die "no macs2 application in path!\n" unless $opts{macs} =~ /\w+/;
	my $qtrack = $self->{qvalue_bdg} || undef;
	return unless $qtrack;
	die "no qvalue bedGraph file $qtrack!\n" unless -e $qtrack;
	my $command = sprintf("%s bdgpeakcall -i %s -c %s -l %s -g %s --no-trackline -o %s ",
		$opts{macs}, $qtrack, $opts{qvalue}, $opts{peaksize}, $opts{peakgap}, 
		$self->{peak}
	);
	my $log = $self->{peak};
	$log =~ s/narrowPeak$/peakcall.out.txt/;
	$command .= " 2> $log";
	return [$command, $self->{peak}, $log];
}

sub generate_bdg2bw_commands {
	my $self = shift;
	my $chromofile = shift;
	my @commands;
	if ($self->{chip_bdg} and $self->{chip_bw}) {
		# we should have both files here
		if (-e $self->{chip_bdg} and -e $self->{chip_bw}) {
			# bigWig already exists so just delete the bdg
			push @commands, [sprintf("rm %s",$self->{chip_bdg}), '', ''];
		}
	}
	if ($self->{lambda_bdg} and $self->{lambda_bw}) {
		if (-e $self->{lambda_bw}) {
			# we must have started with a lambda bigWig so remove the bedGraph
			push @commands, [sprintf("rm %s", $self->{lambda_bdg}), '', ''];
		}
		else {
			my $command = sprintf("%s %s %s %s && rm %s", 
				$opts{wig2bw},
				$self->{lambda_bdg},
				$chromofile,
				$self->{lambda_bw},
				$self->{lambda_bdg},
			);
			push @commands, [$command, $self->{lambda_bw}, '']
		}
	}
	if ($self->{qvalue_bdg} and $self->{qvalue_bw}) {
		die "no wigToBigWig application in path!\n" unless $opts{wig2bw} =~ /\w+/;
		my $command = sprintf("%s %s %s %s && rm %s", 
			$opts{wig2bw},
			$self->{qvalue_bdg},
			$chromofile,
			$self->{qvalue_bw},
			$self->{qvalue_bdg},
		);
		push @commands, [$command, $self->{qvalue_bw}, ''];
	}
	if ($self->{fe_bdg}) {
		# convert this to log2 Fold Enrichment because I like this better
		die "no wigToBigWig application in path!\n" unless $opts{wig2bw} =~ /\w+/;
		die "no manipulate_wig.pl script in path!\n" unless $opts{manwig} =~ /\w+/;
		my $log = $self->{logfe_bw};
		$log =~ s/bw$/out.txt/;
		my $command = sprintf("%s --in %s --log 2 --place 4 --out stdout 2> $log | %s stdin %s %s",
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





