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
use Getopt::Long;
use Parallel::ForkManager;
use Bio::ToolBox 1.67;
use Bio::ToolBox::utility;
use Set::IntervalTree;

our $VERSION = 1.1;

# variables
my $input;
my $outfile;
my $groupfile;
my $cpu;
my $help;

### Documentation
my $docs = <<DOC;

A script to calculate the efficiency of ChIP samples. 

The signal from ChIPs is summed for all identified ChIP peaks, as well as the
remainder of the genome. The fraction of genomic signal corresponding to peaks
is reported. Typically this fraction should be several times higher for ChIPs
than for Input, and this ratio can be interpreted as the efficiency.

Alignment count (point) data is preferred, although fragment data could also be 
used (the numbers will just get very big). Note for count data that the sum of 
on- and off-peak values should be roughly equivalent to the sum of alignments 
used in the analysis, within tolerance of rounding errors for scaled datasets. 

More than one ChIP file can be provided. Each ChIP file is processed in a 
separate thread. 

VERSION: $VERSION

USAGE: get_chip_efficiency.pl -i <peakfile> input.count.bw chip1.count.bw ...

OPTIONS:
    --in <file>             Input file of peak intervals: bed, narrowPeak, etc
    --group <file>          Optional list of groups for each chip file
    --out <file>            The output file, default input basename
    --cpu <int>             Number of threads, default 4
    --help                  Print documentation
DOC

unless (@ARGV) {
	print $docs;
	exit;
}

### Options
GetOptions(
	'in=s'    => \$input,
	'group=s' => \$groupfile,
	'out=s'   => \$outfile,
	'cpu=i'   => \$cpu,
	'help!'   => \$help,
) or die "unrecognized option! See help\n";

# check options
if ($help) {
	print $docs;
	exit;
}

### Check parameters
my @chips = @ARGV;
unless ($input) {
	die " No input peak file provided!\n";
}
unless (@chips) {
	die " One or more ChIP bigWig files must be provided! See help\n";
}
unless ($cpu) {
	$cpu = 4;
}

### Load Input

# Peaks
# we could parse into seqfeatures, but it's not really needed, and this is faster
my $PeakData;
if ( -e $input ) {
	$PeakData = Bio::ToolBox->load_file($input)
		or die " unable to load input file '$input'!\n";
}
else {
	print " no input file '$input'!\n faking it!!!\n";
	$PeakData = Bio::ToolBox->new_data(qw(Chromo Start Stop));
}

# ChIPs
my %chip2group;
if ($groupfile) {
	my $Groups = Bio::ToolBox->load_file($groupfile)
		or die " unable to load group file '$groupfile'!\n";
	$Groups->iterate(
		sub {
			my $row = shift;

			# assuming standard sample and group/condition column
			my $sample = $row->value(0);
			my $group  = $row->value(1);
			$chip2group{$sample} = $group;
		}
	);
}

#### Output
my $Output =
	Bio::ToolBox->new_data(qw(Replicate Dataset OnPeakSum OffPeakSum Efficiency));
$Output->add_comment("Peak file: $input");

### Initialize Fork Manager
my $pm = Parallel::ForkManager->new($cpu);
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) = @_;

		# add child counts to global values
		if ( $exit_code == 0 ) {
			$Output->add_row($data);
		}
		else {
			print " Child $pid failed! See messages\n";
		}
	}
);

### Process peaks
foreach my $chip (@chips) {
	$pm->start and next;
	my $data = collect_peak_data($chip);
	if ( defined $data and $data ) {
		$pm->finish( 0, $data );
	}
	else {
		$pm->finish( 1, undef );
	}
}
$pm->wait_all_children;

### Save
$Output->sort_data( 0, 'i' );    # sort the table by the ChIP names
unless ($outfile) {
	$outfile = $PeakData->path . $PeakData->basename . ".efficiency.txt";
}
my $o = $Output->write_file($outfile);
print " Wrote output file $o\n";

sub collect_peak_data {
	my $file = shift;

	# this is an awkward work around - damned if I do, damned if I don't
	# I need to "verify" the dataset bigwig file if I want to use the simplify dataset
	# name, but the open_db_connection method used by the open_database function
	# doesn't handle "verified" files with a 'file:'prefix. Sigh.
	# this is fixed in Bio::ToolBox 1.68

	my $chipfile = $PeakData->verify_dataset($file);
	unless ($chipfile) {
		print " $file can't be verified!\n";
		return undef;
	}
	my $name  = simplify_dataset_name($chipfile);
	my $group = $chip2group{$name} || $name;

	# Build interval trees
	my %tree;
	$PeakData->iterate(
		sub {
			my $row   = shift;
			my $chr   = $row->seq_id;
			my $start = $row->start - 1;    # use 0-based for this purpose
			my $end   = $row->end;
			$tree{$chr} ||= Set::IntervalTree->new();
			$tree{$chr}->insert( [ $start, $end ], $start, $end );

			# use array reference of start and end for use later
		}
	);

	# open database
	my $db = Bio::ToolBox->open_database($file)
		or warn " Invalid ChIP file specified! Can't open $file!\n";
	return unless ($db);

	# check type of bigWig
	my $type = ref($db);
	unless ( $type eq 'Bio::DB::Big::File' ) {
		print " Currently not supporting $type adaptors! Only Bio::DB::Big bigWig!\n";
		return undef;
	}
	unless ( $db->is_big_wig ) {
		print " database is not a bigWig file!\n";
		return undef;
	}

	# initialize counters
	my $on_peak_sum  = 0;
	my $off_peak_sum = 0;

	# get list of chromosomes
	# this is returned as a hash, so it won't be in same serial order as bigWig file
	# I think that's ok, just maybe not as efficient
	my $chromhash = $db->chroms;

	# walk through each chromosome
	foreach my $chr ( keys %{$chromhash} ) {

		# iterate over the entire length of the chromosome
		# get 500 intervals per iteration to maybe speed things up?
		my $iter =
			$db->get_intervals_iterator( $chr, 0, $chromhash->{$chr}{length}, 500 );
		while ( my $intervals = $iter->next ) {
			foreach my $interval ( @{$intervals} ) {
				my $i_start = $interval->{start};
				my $i_end   = $interval->{end};

				# check for overlapping peaks
				my $results =
					exists $tree{$chr} ? $tree{$chr}->fetch( $i_start, $i_end ) : [];

				if ( scalar @{$results} == 0 ) {

					# wig intervals don't overlap peaks
					# makes things easy
					# multiply value by the length of the interval
					$off_peak_sum += $interval->{value} * ( $i_end - $i_start );
				}
				else {
					# we have overlaps!
					foreach my $r ( @{$results} ) {

						# coordinates for result peak interval
						my $r_start = $r->[0];
						my $r_end   = $r->[1];

						# check base by base
						for ( my $i = $i_start; $i < $i_end; $i++ ) {
							if ( $i >= $r_start and $i < $r_end ) {

								# overlaps a peak
								$on_peak_sum += $interval->{value};
							}
							else {
								# does not overlap a peak
								$off_peak_sum += $interval->{value};
							}
						}
					}
				}
			}
		}
	}

	# reformat to integers, just in case
	if ( $on_peak_sum =~ /\./ ) {
		$on_peak_sum = sprintf( "%.0f", $on_peak_sum );
	}
	if ( $off_peak_sum =~ /\./ ) {
		$off_peak_sum = sprintf( "%.0f", $off_peak_sum );
	}

	# calculate fraction
	my $fraction = sprintf( "%.4f", $on_peak_sum / ( $on_peak_sum + $off_peak_sum ) );

	# finished, return array ref with name, formatted sums, and fraction
	return [
		$name, $group, format_with_commas($on_peak_sum),
		format_with_commas($off_peak_sum), $fraction
	];
}

