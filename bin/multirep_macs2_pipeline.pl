#!/usr/bin/perl

use strict;
use IO::File;
use File::Spec;
use Getopt::Long;

my $VERSION = 1;

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
	fragsize    => 300,
	slocal      => 1000,
	llocal      => 10000,
	qvalue      => 2,
	peaksize    => 300,
	peakgap     => 200,
	targetdep   => 25,
	chrskip     => "chrM|MT|lambda|Adapter|PhiX",
	blacklist   => undef,
	cpu         => 4,
	chipbin     => 10,
	slocalbin   => 50,
	llocalbin   => 100,
	bam2wig     => 'bam2wig.pl',
	bamdedup    => 'bam_partial_dedup.pl',
	macs        => 'macs2',
	manwig      => 'manipulate_wig.pl',
	wig2bw      => 'wigToBigWig',
	bw2bdg      => 'bigWigToBedGraph',
	bedtools    => 'bedtools',
	getdata     => 'get_datasets.pl',
# 	getrdata    => 'get_relative_data.pl',
	printchr    => 'print_chromosome_lengths.pl',
) or die " unrecognized parameter!\n";
$opts{job} = $parallel ? 2 : 1;
my @names;
my @chips;
my @controls;

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
be provided as comma-delimited lists.

Fragment size should be empirically determined by the user, especially when multiple
samples and/or replicates are being used. The same fragment size is used across all
samples and replicates to ensure equal comparisons.

By default, this employs Macs2 local lambda chromatin bias modeling from input. If 
this isn't desired, set slocal and llocal to 0 for both.

Advanced users may provide one processed bigWig file per ChIP or control sample. 


Options:
 Input files
  --chip      file1,file2...    Repeat for each sample set
  --name      text              Repeat for each sample
  --control   file1,file2...    Repeat if matching multiple samples
 
 Output
  --dir       directory         Directory for writing all files ($opts{dir})
  --out       file basename     Base filename for merged output files ($opts{out})
  
 Genome size
  --species   [human,mouse,fish,fly,yeast]   Default ($opts{species})
  --genome    integer           Alternatively give effective genome size
  
 Bam options
  --mapq      integer           Minimum mapping quality, ($opts{species})
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

 Fragment coverage
  --size      integer           Predicted fragment size ($opts{fragsize} bp)
  --slocal    integer           Small local lambda size ($opts{slocal} bp)
  --llocal    integer           Large local lambda size ($opts{llocal} bp)
  --cbin      integer           ChIP fragment bin size ($opts{chipbin} bp)
  --slbin     integer           Small local lambda bin size ($opts{slocalbin} bp)
  --llbin     integer           Large local lambda bin size ($opts{llocalbin} bp)

 Peak calling
  --cutoff    number            Threshold q-value for calling peaks ($opts{qvalue}) 
                                 Higher numbers are more significant, -1*log10(q)
  --tdep      integer           Average sequence depth of bam files in millions ($opts{targetdep})
  --peaksize  integer           Minimum peak size to call ($opts{peaksize} bp)
  --peakgap   integer           Maximum gap between peaks before merging ($opts{peakgap} bp)
  
 Job control
  --cpu       integer           Number of CPUs to use per job ($opts{cpu})
  --job       integer           Number of simultaneous jobs ($opts{job})

 Application  Paths
  --bam2wig   path             ($opts{bam2wig})
  --bamdedup  path             ($opts{bamdedup})
  --macs      path             ($opts{macs})
  --manwig    path             ($opts{manwig})
  --wig2bw    path             ($opts{wig2bw})
  --bedtools  path             ($opts{bedtools})
  --printchr  path             ($opts{printchr})
  --getdata   path             ($opts{getdata})
  
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
	'size=i'                => \$opts{fragsize},
	'slocal=i'              => \$opts{slocal},
	'llocal=i'              => \$opts{llocal},
	'cbin=i'                => \$opts{chipbin},
	'slbin=i'               => \$opts{slocalbin},
	'llbin=i'               => \$opts{llocalbin},
	'cutoff=f'              => \$opts{qvalue},
	'tdep=f'                => \$opts{targetdep},
	'peaksize=i'            => \$opts{peaksize},
	'peakgap=i'             => \$opts{peakgap},
	'cpu=i'                 => \$opts{cpu},
	'job=i'                 => \$opts{job},
	'bam2wig=s'             => \$opts{bam2wig},
	'bamdedup=s'            => \$opts{bamdedup},
	'macs=s'                => \$opts{macs},
	'manwig=s'              => \$opts{manwig},
	'wig2bw=s'              => \$opts{wig2bw},
	'bedtools=s'            => \$opts{bedtools},
	'getdata=s'             => \$opts{getdata},
	'printchr=s'            => \$opts{printchr},
) or die "unrecognized option(s)!\n";



### Begin main pipeline
my $start = time;
my @finished_commands;
check_inputs();
print_start();
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
}

sub print_start {
	print "======== ChIPSeq multi-replicate pipeline ==========\n";
	print "\nversion $VERSION\n";
	print "\n\n======= Samples\n";
	for my $i (0 .. $#names) {
		printf " %s: %s\n", $names[$i], $chips[$i]; 
		printf " Control: %s\n", $controls[$i] if $i < scalar(@controls);
	}
	print "\n\n======= Configuration\n";
	foreach my $k (sort {$a cmp $b} keys %opts) {
		printf "%10s %30s\n", $k, $opts{$k};
	}
	print "\n\n";
}

sub generate_job_file_structure {
	# we pass this off to ChIPjob package with 3 options: name, chip files, control files
	my @jobs;
	# check for a unversal control
	if (scalar(@controls) == 1 and scalar(@chips) > 1) {
		print "Using only one control for multiple ChIP experiments\n";
		my $universal_control = shift @controls;
		foreach (@names) {
			push @controls, 'universal';
		}
		push @jobs, ChIPjob->new($opts{out} . "_control", '', $universal_control);
	}
	for my $i (0 .. $#names) {
		push @jobs, ChIPjob->new($names[$i], $chips[$i], $controls[$i] || undef);
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
			print "=== Job: $command\n";
			$pm->start and next;
			# in child
			system($command);
			$pm->finish;
		}
		$pm->wait_all_children;
	}
	else {
		foreach my $command (@$commands) {
			print "=== Job: $command\n\n";
			system($command);
		}
	}
	push @finished_commands, @$commands;
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
	if ($opts{slocal} and $opts{llocal}){
		my @commands;
		foreach my $Job (@Jobs) {
			push @commands, $Job->generate_lambda_control_commands;
		}
		if (@commands) {
			print "\n\n======= Generate control lambda files\n";
			execute_commands(\@commands);
		}
	}
	else {
		# we don't want lambda control, so just rename the file
		my @commands;
		foreach my $Job (@Jobs) {
			if ($Job->{d_control_bdg}) {
				push @commands, sprintf("mv %s %s", $Job->{d_control_bdg}, $Job->{lambda_bdg});
				$Job->{lambda_bdg} =~ s/\.lambda_control// if $Job->{lambda_bdg};
				$Job->{lambda_bw} =~ s/\.lambda_control// if $Job->{lambda_bw};
			}
		}
		if (@commands) {
			print "\n\n======= Skipping control lambda\n";
			execute_commands(\@commands);
		}
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
	my $chromofile = generate_chr_file();
	my @commands;
	foreach my $Job (@Jobs) {
		push @commands, $Job->generate_bdg2bw_commands($chromofile);
	}
	print "\n\n======= Converting bedGraph files\n";
	execute_commands(\@commands);
	unlink $chromofile;
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
	print "\n\n======= Merging called narrowPeak files\n";
	my $command = "cat ";
	foreach my $Job (@Jobs) {
		if ($Job->{peak}) {
			$command .= sprintf("%s ", $Job->{peak});
		}
	}
	$command .= sprintf(" | %s sort -i - | %s merge -i - > %s",
		$opts{bedtools}, $opts{bedtools}, 
		File::Spec->catfile($opts{dir}, $opts{out} . '.bed') );
	execute_commands([$command]);
}

sub run_rescore {
	print "\n\n======= Re-scoring all merged peaks\n";
	
	my $input = File::Spec->catfile($opts{dir}, $opts{out} . '.bed');
	die "unable to find merged bed file '$input'!\n" unless -e $input;
	my $output1 = File::Spec->catfile($opts{dir}, $opts{out} . '_qvalue.txt');
	my $output2 = File::Spec->catfile($opts{dir}, $opts{out} . '_log2FE.txt');
	
	my $command1 = sprintf("%s --method mean --cpu %s --in %s --out %s ",
		$opts{getdata}, $opts{cpu}, $input, $output1);
	my $command2 = sprintf("%s --method mean --cpu %s --in %s --out %s ",
		$opts{getdata}, $opts{cpu}, $input, $output2);
	foreach my $Job (@Jobs) {
		if ($Job->{qvalue_bw}) {
			$command1 .= sprintf("--data %s ", $Job->{qvalue_bw});
		}
		if ($Job->{logfe_bw}) {
			$command2 .= sprintf("--data %s ", $Job->{logfe_bw});
		}
	}
	$command1 .= sprintf(" 2>&1 > %s", File::Spec->catfile($opts{dir}, 
		'qvalue_get_datasets.out.txt'));
	$command2 .= sprintf(" 2>&1 > %s", File::Spec->catfile($opts{dir}, 
		'log2FE_get_datasets.out.txt'));
	execute_commands([$command1, $command2]);
}

sub finish {
	my @combined_output;
	foreach my $c (@finished_commands) {
		if ($c =~ m/> ([^ ]+\.out\.txt)/) {
			my $f = $1;
			push @combined_output, "===Job: $c\n";
			unless (-z $f) {
				my $fh = IO::File->new($f) or next;
				push @combined_output, <$fh>;
				push @combined_output, "\n\n";
				$fh->close;
			}
			unlink $f if -e _;
		}
	}
	my $file = File::Spec->catfile($opts{dir}, $opts{out} . "_job_output_logs.txt");
	my $fh = IO::File->new($file, "w");
	foreach (@combined_output) {
		$fh->print($_);
	}
	$fh->close;
	print "\n Combined all job output log files into '$file'\n";
	printf " Finished in %.1f minutes\n", (time -$start) / 60;
	print "\n======== Finished ChIPSeq multi-replicate pipeline ==========\n";
}



############### ChIPJob Package ########################################################

package ChIPjob;

sub new {
	# pass name, comma-list of ChIP files, and comma-list of control files
	my ($class, $name, $chip, $control) = @_;
	my $namepath = File::Spec->catfile($opts{dir}, $name);
	my $self = {
	    name => $name,
	    chip_bams => undef,
	    control_bams => undef,
	    chip_dedup_bams => [],
	    chip_use_bams => [],
	    control_dedup_bams => [],
	    control_use_bams => [],
	    chip_bw => undef,
	    chip_bdg => undef,
	    d_control_bdg => undef,
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
	};
	
	# check the ChIP files
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
	if ($chip =~ /\.bam$/i) {
		$self->{chip_bams} = [split(',', $chip)];
		$self->{chip_bdg} = "$namepath.fragment.bdg";
		$self->{chip_bw} = "$namepath.fragment.bw";
		$self->{qvalue_bdg} = "$namepath.qvalue.bdg";
		$self->{qvalue_bw} = "$namepath.qvalue.bw";
		$self->{fe_bdg} = "$namepath.FE.bdg";
		$self->{logfe_bw} = "$namepath.log2FE.bw";
		$self->{peak} = "$namepath.narrowPeak";
	}
	else {
		# must be a control
		$self->{chip_bams} = [];
	}
	
	# check the control files
	if ($control eq 'universal') {
		# this is just a ChIP Job using a common universal control lambda
		$self->{lambda_bdg} = File::Spec->catfile($opts{dir}, 
			$opts{out} . "_control.lambda_control.bdg");
	}
	elsif ($control =~ /\.(?:bw|bigwig)$/i) {
		die "only one control biwWig file is allowed per experiment!\n" if $chip =~ /,/;
		$self->{lambda_bw} = $control;
		$self->{lambda_bdg} = $control;
		$self->{lambda_bdg} =~ s/(?:bw|bigwig)$/bdg/i;
	}
	elsif ($control =~ /\.bam$/i) {
		# we have control bams to process
		$self->{control_bams} = [split(',', $control)];
		$self->{lambda_bdg} = "$namepath.lambda_control.bdg";
		$self->{lambda_bw} = "$namepath.lambda_control.bw";
		$self->{d_control_bdg} = "$namepath.dlocal.bdg";
		$self->{s_control_bdg} = "$namepath.slocal.bdg";
		$self->{l_control_bdg} = "$namepath.llocal.bdg";
		$self->{sl_control_bdg} = "$namepath.sllocal.bdg";
		$self->{sld_control_bdg} = "$namepath.sldlocal.bdg";
	}
	else {
		# must be just a chip without corresponding control
		$self->{control_bams} = [];
	}
	return bless $self, $class;
}

sub generate_dedup_commands {
	my $self = shift;
	my @commands;
	if (defined $self->{chip_bams}) {
		for (my $i = 0; $i < scalar @{$self->{chip_bams}}; $i++) {
			my $in = $self->{chip_bams}->[$i];
			my (undef,undef,$out) = File::Spec->splitpath($in);
			$out = File::Spec->catfile($opts{dir}, $out);
			$out =~ s/\.bam$/.dedup.bam/i;
			$self->{chip_dedup_bams}->[$i] = $out;
			my $command = sprintf "%s --in %s --out %s ", 
				$opts{bamdedup},
				$in,
				$out;
			if ($opts{dupfrac} > 0) {
				$command .= sprintf("--random --seed 1 --frac %s ", $opts{dupfrac});
			}
			if (defined $opts{maxdup}) {
				$command .= sprintf("--max %s ", $opts{maxdup});
			}
			if ($opts{paired}) {
				$command .= "--pe ";
			}
			my $log = $out;
			$log =~ s/\.bam$/.out.txt/i;
			$command .= " 2>&1 > $log";
			push @commands, $command;
		}
	}
	if (defined $self->{control_bams}) {
		for (my $i = 0; $i < scalar @{$self->{control_bams}}; $i++) {
			my $in = $self->{control_bams}->[$i];
			my (undef,undef,$out) = File::Spec->splitpath($in);
			$out = File::Spec->catfile($opts{dir}, $out);
			$out =~ s/\.bam$/.dedup.bam/i;
			$self->{control_dedup_bams}->[$i] = $out;
			my $command = sprintf "%s --in %s --out %s ", 
				$opts{bamdedup},
				$in,
				$out;
			if ($opts{dupfrac} > 0) {
				$command .= sprintf("--random --seed 1 --frac %s ", $opts{dupfrac});
			}
			if (defined $opts{maxdup}) {
				$command .= sprintf("--max %s ", $opts{maxdup});
			}
			if ($opts{paired}) {
				$command .= "--pe ";
			}
			my $log = $out;
			$log =~ s/\.bam$/.out.txt/i;
			$command .= " 2>&1 > $log";
			push @commands, $command;
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
	my @commands;
	
	# ChIP bams
	if (scalar @{$self->{chip_use_bams}}) {
		# we have bam files to convert to bw
		# base command
		my $command = sprintf(
			"%s --in %s --out %s --qual %s --bin %s --cpu %s --rpm --mean --bw --bwapp %s ", 
			$opts{bam2wig}, 
			join(',', @{$self->{chip_use_bams}}), 
			$self->{chip_bw},
			$opts{mapq},
			$opts{chipbin},
			$opts{cpu},
			$opts{wig2bw}
		);
		# paired or single options
		if ($opts{paired}) {
			$command .= sprintf("--pe --span --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
		}
		else {
			$command .= sprintf("--extend --extval %s ", $opts{fragsize});
		}
		# additional filtering
		if ($opts{blacklist}) {
			$command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		my $log = $self->{chip_bw};
		$log =~ s/bw$/out.txt/;
		$command .= " 2>&1 > $log";
		# done
		push @commands, $command;
	}
	
	# Control bams
	if (scalar @{$self->{control_use_bams}}) {
		# base command
		my $command = sprintf(
			"%s --in %s --qual %s --cpu %s --rpm --mean --cspan --bdg ", 
			$opts{bam2wig}, 
			join(',', @{$self->{control_use_bams}}), 
			$opts{mapq},
			$opts{cpu}
		);
		# paired options
		if ($opts{paired}) {
			$command .= sprintf("--pe --minsize %s --maxsize %s ", 
				$opts{minsize}, $opts{maxsize});
		}
		# additional filters
		if ($opts{blacklist}) {
			$command .= sprintf("--blacklist %s ", $opts{blacklist});
		}
		if ($opts{chrskip}) {
			$command .= sprintf("--chrskip \'%s\' ", $opts{chrskip});
		}
		
		# now duplicate the base command for each size lambda control
		
		# d control
		my $command1 = $command;
		$command1 .= sprintf("--extval %s --bin %s --out %s ", 
			$opts{fragsize}, $opts{chipbin}, $self->{d_control_bdg}
		);
		my $log = $self->{d_control_bdg};
		$log =~ s/bdg$/out.txt/;
		$command1 .= " 2>&1 > $log";
		push @commands, $command1;
		
		# small local lambda
		if ($opts{slocal}) {
			my $command2 = $command;
			$command2 .= sprintf("--extval %s --scale %.6f --bin %s --out %s ",
				$opts{slocal}, $opts{fragsize}/$opts{slocal}, $opts{slocalbin}, 
				$self->{s_control_bdg});
			$log = $self->{s_control_bdg};
			$log =~ s/bdg$/out.txt/;
			$command2 .= " 2>&1 > $log";
			push @commands, $command2;
		}
		# large local lambda
		if ($opts{llocal}) {
			my $command3 = $command;
			$command3 .= sprintf("--extval %s --scale %.6f --bin %s --out %s ",
				$opts{llocal}, $opts{fragsize}/$opts{llocal}, $opts{llocalbin}, 
				$self->{l_control_bdg});
			$log = $self->{l_control_bdg};
			$log =~ s/bdg$/out.txt/;
			$command3 .= " 2>&1 > $log";
			push @commands, $command3
		}
	}
	
	# finished
	return @commands;
}

sub generate_lambda_control_commands {
	my $self = shift;
	my $dfile = $self->{d_control_bdg};
	my $sfile = $self->{s_control_bdg};
	my $lfile = $self->{l_control_bdg};
	return unless ($dfile and $sfile and $lfile);
	die "no d control bedGraph file '$dfile'!\n" unless (-e $dfile);
	die "no small control bedGraph file '$sfile'!\n" unless (-e $sfile);
	die "no large control bedGraph file '$lfile'!\n" unless (-e $lfile);
	
	my $log = $self->{lambda_bdg};
	$log =~ s/bdg$/out.txt/;
	
	# first step
	my $command = sprintf("%s bdgcmp -m max -t $sfile -c $lfile -o %s ", 
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
	
	return $command;
}

sub convert_bw_to_bdg {
	my $self = shift;
	my @commands;
	if ($self->{chip_bw} and -e $self->{chip_bw}) {
		push @commands, sprintf("%s %s %s", $opts{bw2bdg}, $self->{chip_bw}, 
			$self->{chip_bdg});
	}
	if ($self->{lambda_bw} and -e $self->{lambda_bw}) {
		push @commands, sprintf("%s %s %s", $opts{bw2bdg}, $self->{lambda_bw}, 
			$self->{lambda_bdg});
	}
	return @commands;
}

sub generate_enrichment_commands {
	my $self = shift;
	my $chip = $self->{chip_bdg} || undef;
	my $lambda = $self->{lambda_bdg} || undef;
	return unless ($chip and $lambda);
	die "no ChIP fragment file $chip!\n" unless -e $chip;
	die "no control lambda fragment file $lambda!\n" unless -e $lambda;
	
	my $command = sprintf("%s bdgcmp -t %s -c %s -S %s -m qpois FE -o %s %s ", 
		$opts{macs}, $chip, $lambda, $opts{targetdep}, $self->{qvalue_bdg}, 
		$self->{fe_bdg});
	my $log = $self->{qvalue_bdg};
	$log =~ s/bdg$/out.txt/;
	$command .= " 2> $log ";
	return $command;
}

sub generate_peakcall_commands {
	my $self = shift;
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
	return $command;
}

sub generate_bdg2bw_commands {
	my $self = shift;
	my $chromofile = shift;
	my @commands;
	if ($self->{chip_bdg} and $self->{chip_bw}) {
		# we should have both files here
		if (-e $self->{chip_bdg} and -e $self->{chip_bw}) {
			# bigWig already exists so just delete the bdg
			push @commands, sprintf("rm %s",$self->{chip_bdg});
		}
	}
	if ($self->{lambda_bdg} and $self->{lambda_bw}) {
		if (-e $self->{lambda_bw}) {
			# we must have started with a lambda bigWig so remove the bedGraph
			push @commands, sprintf("rm %s", $self->{lambda_bdg});
		}
		else {
			push @commands, sprintf("%s %s %s %s && rm %s", 
				$opts{wig2bw},
				$self->{lambda_bdg},
				$chromofile,
				$self->{lambda_bw},
				$self->{lambda_bdg},
			);
		}
	}
	if ($self->{qvalue_bdg} and $self->{qvalue_bw}) {
		my $command = sprintf("%s %s %s %s && rm %s", 
			$opts{wig2bw},
			$self->{qvalue_bdg},
			$chromofile,
			$self->{qvalue_bw},
			$self->{qvalue_bdg},
		);
		push @commands, $command;
	}
	if ($self->{fe_bdg}) {
		# convert this to log2 Fold Enrichment because I like this better
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
		push @commands, $command;
	}
	return @commands;
}





