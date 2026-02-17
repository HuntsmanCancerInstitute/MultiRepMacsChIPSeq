#!/usr/bin/perl

# Timothy J. Parnell, PhD
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112
#
# This package is free software; you can redistribute it and/or modify
# it under the terms of the Artistic License 2.0.
#

use warnings;
use strict;
use English qw(-no_match_vars);
use Getopt::Long;
use File::Copy;
use File::Which;
use List::Util qw(sum0);
use Bio::ToolBox 2.03;
use Bio::ToolBox::utility qw(format_with_commas);

our $VERSION = 0.1;

# user variables
my $rmsk_file;
my $in_file;
my $out_file;
my @datasets;
my $min_length = 1000;
my $min_zscore = 1;
my $chrskip    = 'chrM|MT|alt|chrUn|random|EBV|Adapter|Lambda|PhiX';
my $save       = 0;
my $cpu        = 4;
my $getdata    = which('get_datasets.pl');
my $mandata    = which('manipulate_datasets.pl');
my $help;

### Documentation
my $docs = <<DOC;

A script to generate a new exclusion interval file based on extreme coverage
over known repetitive elements.

Repetitive elements are a known potential source of false-positive peak calls.
However, not all repetitive elements result in peaks. Most elements are
degenerate enough or have enough copies scattered across the genome that
they can be ignored by aligners and not generate significant coverage to be
a problem. But a few can generate extraordinarily high coverage that can
be called as peaks, detracting signal from real peaks. Often these are long,
complete transposable elements or centromeric elements.

While publicly available exclusion lists are available (so called "black"
lists), differences in lab strains, conditions, and aligners only make these
partially useful. Ideally, exclusions would be generated from Input chromatin
coverage, but this is not always available, particularly for ATAC-Seq,
Cut&Run, Cut&Tag, and related experiments.

This script takes a public list of known repetitive elements for a genome
(RepeatMasker files obtained from UCSC), filters for single or adjoining
elements that exceed the minimum length, scores them for normalized counts,
and identifies those with extreme coverage to be used as exclusion
intervals for bam filtering. Typically in empirical testing, the
distribution of coverage is fairly narrow, with only a tiny fraction of
elements having extreme coverage, which means having a relatively low
threshold (1 standard deviation from the mean) is sufficient for calling
outliers.

In general, all sample and replicate bam files in an experiment should be
provided. The average is taken across the samples to ensure that only those
elements with consistently high coverage are used for exclusion.

VERSION: $VERSION

USAGE: rmsk2exclusion.pl -r rmsk.txt -o exclusion.bed  File1.bam  File2.bam ...

OPTIONS:

Input:
  -r --rmsk <file>         Input RepeatMasker file. Required. These can often
                            be obtained for your genome from
                            https://genome.ucsc.edu
  -i --input <file>        Bed file of custom repetitive elements to score.
                            This can be used as an alternative --rmsk.
  -o --out <file>          Output file basename. Required.
  -d --data <file>         Specify one or more bam files to score. This option
                            may be repeated for each bam file, or bam files
                            may be simply appended to the end of the command.

Filtering Options:
  -l --length <int>        Minimum length of repeat element to consider ($min_length bp)
  -z --zscore <float>      Minimum z-score for scoring elements as exclusion ($min_zscore)
  --chrskip <regex>        Chromosome skip regex ($chrskip)
  

General:
 --save                    Save the intermediate scored repeat element file.
 -c --cpu <int>            Number of threads to use
 --getdata <path>          Path to get_datasets.pl ($getdata)
 --mandata <path>          Path to manipulate_datasets.pl ($mandata)
 -h --help                 This help



DOC

# Options
unless (@ARGV) {
	print $docs;
	exit;
}
GetOptions(
	'r|rmsk=s'   => \$rmsk_file,
	'i|in=s'     => \$in_file,
	'o|out=s'    => \$out_file,
	'd|data=s'   => \@datasets,
	'l|length=i' => \$min_length,
	'z|zscore=f' => \$min_zscore,
	'chrskip=s'  => \$chrskip,
	'save!'      => \$save,
	'c|cpu=i'    => \$cpu,
	'getdata=s'  => \$getdata,
	'mandata=s'  => \$mandata,
	'h|help!'    => \$help,
) or die "unrecognized option!\n";

if ($help) {
	print $docs;
	exit;
}

# Check input options
unless ( $rmsk_file or $in_file ) {
	print "FATAL: No input file!\n";
	exit 1;
}
unless ($out_file) {
	print "FATAL: No output file name provided!\n";
	exit 1;
}
unless (@datasets) {
	if (@ARGV) {
		@datasets = @ARGV;
	}
	else {
		print "FATAL: No dataset bam files provided!\n";
		exit 1;
	}
}
unless ($getdata) {
	print "FATAL: no get_datasets.pl application path provided or found!\n";
	exit 1;
}
unless ($mandata) {
	print "FATAL: no manipulate_datasets.pl application path provided or found!\n";
	exit 1;
}

# temporary output files
my $rmsk_bed_file = sprintf( "rmsk.%s.bed",          $PID );
my $rmsk_scr_file = sprintf( "rmsk_score.%s.txt.gz", $PID );

# global variables
my @zcols;

#### main
if ($in_file) {
	load_input_file();
}
elsif ($rmsk_file) {
	load_rmsk_file();
}
score_repeats();
filter_repeats();

# finish
unlink $rmsk_bed_file;
if ($save) {
	$out_file =~ s/ \.bed (?:\.gz)? $//xi;
	$out_file .= '.scored.txt.gz';
	copy( $rmsk_scr_file, $out_file );
	printf " Saved intermediate scored interval file '%s'\n", $out_file;
}
else {
	unlink $rmsk_scr_file;
}
exit;

#### Functions
sub load_input_file {

	# open the file
	printf " Loading %s...\n", $in_file;
	my $Stream = Bio::ToolBox->load_file( in => $in_file, stream => 1 )
		or die "unable to open file '$in_file'!";
	unless ( $Stream->feature_type eq 'coordinate' ) {
		print "Input file '$in_file' does not appear to have coordinates\n";
		exit 1;
	}

	# convert to a simple sorted bed
	my $Data;
	my $name_i = $Stream->name_column;
	if ($name_i) {
		$Data = Bio::ToolBox->new_bed(4);
	}
	else {
		$Data = Bio::ToolBox->new_bed(3);
	}
	while ( my $row = $Stream->next_row ) {
		next if $row->seq_id =~ /$chrskip/xi;
		my @bits = ( $row->seq_id, $row->start - 1, $row->end );
		if ($name_i) {
			push @bits, $row->value($name_i);
		}
		$Data->add_row( \@bits );
	}
	$Stream->close_fh;
	unless ( $Data->number_rows ) {
		print " ERROR: No intervals loaded!\n";
		exit 1;
	}

	# save
	$Data->gsort_data;
	printf " Using %s elements\n", format_with_commas( $Data->number_rows );
	$Data->write_file($rmsk_bed_file)
		or die "unable to write temporary bed file $rmsk_bed_file!";
}

sub load_rmsk_file {

	# open the file as a stream to avoid excessive memory usage
	my $Stream = Bio::ToolBox->load_file( in => $rmsk_file, stream => 1 )
		or die "Unable to open file '$rmsk_file'!";
	printf " Loading %s...\n", $rmsk_file;

	# find columns
	my $chr_i   = $Stream->find_column('genoName|chromosome');
	my $start_i = $Stream->find_column('genoStart|start');
	my $stop_i  = $Stream->find_column('genoEnd|end|stop');
	my $name_i  = $Stream->find_column('repClass|class|family');

	unless ( $chr_i and $start_i and $stop_i ) {
		printf
"FATAL: Unable to identify genome chromosome, start, and stop columns in %s!\n",
			$rmsk_file;
		exit 1;
	}

	# generate simple bed first
	# we need to sort the elements and can't do that until we have standard coordinates
	# RepeatMasker files from UCSC are rarely if ever sorted in a sane fashion
	my $Data = Bio::ToolBox->new_bed(4);
	while ( my $row = $Stream->next_row ) {
		next if $row->value($chr_i) =~ /$chrskip/xi;
		$Data->add_row(
			[
				$row->value($chr_i),
				$row->value($start_i),
				$row->value($stop_i),
				$row->value($name_i)
			]
		);
	}
	$Stream->close_fh;
	unless ( $Data->number_rows ) {
		print " ERROR: No features loaded!\n";
		exit 1;
	}
	printf " Loaded %s features from %s\n", format_with_commas( $Data->number_rows ),
		$rmsk_file;

	# seed the first interval
	printf " Sorting...\n";
	$Data->gsort_data;
	my $current = $Data->get_row(1)->row_values;

	# pull adjoining intervals that exceed minimum length
	my $Final = Bio::ToolBox->new_bed(4);
	$Data->iterate(
		sub {
			my $row = shift;
			if ( $row->seq_id eq $current->[0] and $row->start <= $current->[2] ) {

				# update current
				$current->[2] = $row->end - 1;
				$current->[3] .= sprintf ";%s", $row->name;
			}
			else {
				if ( ( $current->[2] - $current->[1] ) >= $min_length ) {
					$Final->add_row($current);
				}
				$current = $row->row_values;
			}
		}
	);

	# write last one
	if ( ( $current->[2] - $current->[1] ) >= $min_length ) {
		$Final->add_row($current);
	}

	# finish
	printf " Selected %s merged features to score\n",
		format_with_commas( $Final->number_rows );
	$Final->write_file($rmsk_bed_file)
		or die "unable to write temporary bed file $rmsk_bed_file!";
}

sub score_repeats {

	# collect read counts as FPKM value over each element
	# normalize to all counts in the genome
	my $command1 = sprintf
		"%s --in %s --noparse --method count --fpkm genome --cpu %s --gz --out %s %s",
		$getdata,
		$in_file || $rmsk_bed_file,
		$cpu,
		$rmsk_scr_file,
		join( q( ), map {"--data $_"} @datasets );

	printf " Scoring features with %d bam files...\n", scalar(@datasets);
	system($command1) == 0
		or die " FATAL: execution '$command1' failed!\n";

	# check the columns of the output file
	my $Stream = Bio::ToolBox->load_file( in => $rmsk_scr_file, stream => 1 );
	my @ids;
	for my $i ( 1 .. $Stream->number_columns ) {
		if ( $Stream->name($i) =~ /_FPKM$/ ) {
			push @ids, $i;
		}
	}
	$Stream->close_fh;

	# calculate the Z-score for each provided replicate, record as new column
	my $command2 = sprintf "%s --in %s --func zscore --place n --index %s",
		$mandata, $rmsk_scr_file, join( ',', @ids );
	print "\n Calculating Z-scores...\n";
	system($command2) == 0
		or die " FATAL: execution '$command2' failed!\n";

	# remember the new Z-score columns
	@zcols = ( $ids[-1] + 1 .. $ids[-1] + scalar(@ids) );
}

sub filter_repeats {

	# load the scored repeat element file
	printf "\n Filtering for outlier features...\n";
	my $Stream = Bio::ToolBox->load_file( in => $rmsk_scr_file, stream => 1 );

	# prepare final selected Data file
	my $Final = Bio::ToolBox->new_bed(5);

	# select
	my $n     = scalar(@zcols);
	my $count = 0;
	while ( my $row = $Stream->next_row ) {
		my $score = sum0( map { $row->value($_) } @zcols ) / $n;
		if ( $score >= $min_zscore ) {
			$count++;
			my @bits = (
				$row->seq_id, $row->start - 1, $row->end,
				$row->name || sprintf( "region%d", $count ), sprintf( "%.3f", $score )
			);
			$Final->add_row( \@bits );
		}
	}
	$Stream->close_fh;

	# finish
	printf " Retained %s intervals as exclusion regions\n", format_with_commas($count);
	my $success = $Final->write_file($out_file);
	printf " Wrote final exclusion file '%s'\n", $success;
}

