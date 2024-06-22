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
use Bio::ToolBox 1.70;
use List::Util qw(max sum0);
use Statistics::Lite qw(median);

our $VERSION = 1.2;

unless (@ARGV) {
	print <<USAGE;

A script to combine replicate values into a single value and 
written as a new data file. Replicates may be combined with a 
number of different methods. Replicate groups are specified with 
a sample file: a two-column text file with sample identifiers in 
the first column and its group identifier in the second. Sample 
identifiers must match the column name in the input replicate file. 
File compression is natively handled.

Version: $VERSION

Usage: combine_replicate_data.pl -i rep_counts.txt -s samples.txt -o sample_counts.txt

Required:
  -i --in <file>          Input file of replicate counts
  -s --sample <file>      File of replicate samples and groups
  -o --out <file>         Output file name

Optional:
  -m --method <text>      Method of combining: sum mean median max (mean)
  -f --format <integer>   Number of decimals to format combined values (0)

USAGE
	exit;
}

### Options
my $infile;
my $samplefile;
my $method = 'mean';
my $format = 0;
my $outfile;
GetOptions(
	'i|in=s'     => \$infile,
	's|sample=s' => \$samplefile,
	'm|method=s' => \$method,
	'f|format=i' => \$format,
	'o|out=s'    => \$outfile
) or die "unrecognized option!\n";

die "no input file specified!\n"  unless $infile;
die "no sample file specified!\n" unless $samplefile;
die "no output file specified!\n" unless $outfile;

# methods
my $method_sub;
if ( $method eq 'mean' ) {
	$method_sub = sub { return sum0(@_) / scalar(@_) };
}
elsif ( $method eq 'sum' ) {
	$method_sub = \&sum0;
}
elsif ( $method eq 'max' ) {
	$method_sub = \&max;
}
elsif ( $method eq 'median' ) {
	$method_sub = \&median;
}

# sprintf string for formatting
my $printer = '%.' . $format . 'f';

my $start_time = time;

### Input file stream
my $InData = Bio::ToolBox->load_file(
	file   => $infile,
	stream => 1,
) or die "unable to open '$infile'!\n";

### Sample file
my $Sample = Bio::ToolBox->load_file($samplefile)
	or die "unable to open '$samplefile'!\n";
unless ( $Sample->number_columns == 2 ) {
	die "sample file has more than two columns!\n";
}

### Identify groups
my %group2rep;
my $repnumber = 0;
$Sample->iterate(
	sub {
		my $row = shift;

		# assuming sample group are the two columns
		# find the associating column in the input Data for each replicate
		my $sample = $row->value(1);
		my $group  = $row->value(2);
		my $i      = $InData->find_column("^$sample\$");
		unless ( defined $i ) {
			warn " can't find column for sample '$sample'!\n";
			next;
		}
		$group2rep{$group} ||= [];
		push @{ $group2rep{$group} }, $i;
		$repnumber++;
	}
);

# put into array of arrays
my @names  = sort { $a cmp $b } keys %group2rep;
my @groups = map  { $group2rep{$_} } @names;

### Output data stream
# first identify which columns we've got
my @standard;
{    # naked block
	my $chrom_i = $InData->chromo_column;
	my $start_i = $InData->start_column;
	my $stop_i  = $InData->stop_column;
	my $strnd_i = $InData->strand_column;
	my $id_i    = $InData->id_column;
	my $name_i  = $InData->name_column;
	my $type_i  = $InData->type_column;
	push @standard, $chrom_i if ( defined $chrom_i );
	push @standard, $start_i if ( defined $start_i );
	push @standard, $stop_i  if ( defined $stop_i );
	push @standard, $strnd_i if ( defined $strnd_i );
	push @standard, $id_i    if ( defined $id_i );
	push @standard, $name_i  if ( defined $name_i );
	push @standard, $type_i  if ( defined $type_i );
}

# assemble the final array of indices for the final file
# sort the columns, these should already be in order, but just in case
@standard = sort { $a <=> $b } @standard;

my @standard_names = map { $InData->name($_) } @standard;
my $OutData        = Bio::ToolBox->new_data(
	columns => [ @standard_names, @names ],
	stream  => 1,
	out     => $outfile,
);
$OutData->add_comment("$method sample data from $infile");

### Combine data and write to output
while ( my $row = $InData->next_row ) {

	# first collect the standard data
	my @data = map { $row->value($_) } @standard;

	# then collect
	foreach my $g (@groups) {
		push @data, sprintf( $printer, &{$method_sub}( map { $row->value($_) } @{$g} ) );
	}

	# store data
	$OutData->add_row( \@data );
}

### Finish
$InData->close_fh;
$OutData->close_fh;
printf(
	" Combined %d replicates in %d groups\n Wrote file $outfile in %.1f seconds\n",
	$repnumber, scalar(@groups), time - $start_time
);

