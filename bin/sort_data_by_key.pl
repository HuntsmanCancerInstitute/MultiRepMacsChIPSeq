#!/usr/bin/env perl

# documentation at end of file

use warnings;
use strict;
use English qw(-no_match_vars);
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use List::Util qw(sum0);
use Scalar::Util qw(looks_like_number);
use Bio::ToolBox::Data;
use Bio::ToolBox::utility qw(parse_list);

our $VERSION = 0.2;

my $docs = <<DOC;

This script will sort the features in a data file by the groups
indicated in a second provided file. Features are sorted numerically
by the row mean values within each group. The primary intention
is to improve plotting of data by extrinsic grouping rather than
intrinsic ordering.

The data file should be a table, such as generated by `get_datasets.pl`
or `get_relative_data.pl`, where the ID and/or Name column is a unique
identifier and multiple data columns of numeric data exist. 

The key or matrix file should contain the same unique feature identifiers
(Name or ID) in the first column as the input data file. The groups should
be in one or more subsequent columns. If more than one column, the values 
are concatenated with a '-' delimiter. The boolean intersection matrix 
file from `intersect_peaks.pl` may be used here, for example. 

An opportunity is provided to normalize the data values by specifying 
the reference columns. In each row, the mean value from the reference 
columns is subtracted from all of the data values. The purpose is to
generate a relative difference table as a convenience.

Multiple output files are written:

	- A sorted data table with a basename suffix of '.sorted'
	- A '.row_groups.txt' annotation file with the row groups
	- A summary data file for each group for plotting mean profiles

VERSION: $VERSION

USAGE: sort_data_by_key.pl  -i <dataset_file>  -k <matrix_file> 

REQUIRED:
  -i --input <file>      The input data file
  -k --key <file>        The matrix file of keys and groups

OPTIONS:
  -g --group <index>     The group column(s) in the key file
                           default is all columns except the first ID column
  -v --values <index>    Specific data column indexes to sort by
                           default is all identifiable data columns
  -n --norm <index>      Specify the reference data columns for normalization
  -f --format <int>      Specify the decimal positions when normalizing
  -s --sort [count|name] Specify how to sort the groups: by decreasing 
                           count or by name. Default is count.
  -h --help              Show this help                           

DOC

my $data_file;
my $matrix_file;
my $grouplist;
my $valuelist;
my $normlist;
my $format;
my $sort = 'count';
my $help;


unless (@ARGV) {
	print $docs;
	exit;
}
GetOptions(
	'i|input=s'  => \$data_file,
	'k|key=s'    => \$matrix_file,
	'g|group=s'  => \$grouplist,
	'v|value=s'  => \$valuelist,
	'n|norm=s'   => \$normlist,
	'f|format=i' => \$format,
	's|sort=s'   => \$sort,
	'h|help!'    => \$help,
) or die "unrecognized option!\n";

# check options
if ($help) {
	print $docs;
	exit;
}
unless ($data_file) {
	print " ERROR! No input data file specified!\n";
	exit 1;
}
unless ($matrix_file) {
	print " ERROR! No key matrix file specified!\n";
	exit 1;
}
if (defined $format) {
	$format = '%.' . $format . 'f';
}


# global hashes
my %name2group;
my %name2row;
my %group2count;
my %group2data;

# normalization indexes for subtracting
my @norm_indexes;
if ($normlist) {
	@norm_indexes = parse_list($normlist);
}

# group index
my @group_indexes;
if ($grouplist) {
	@group_indexes = parse_list($grouplist);
}

# value index
my @val_indexes;
if ($valuelist) {
	@val_indexes = parse_list($valuelist);
}

# Load matrix
print " Loading matrix $matrix_file...\n";
my $Matrix = Bio::ToolBox::Data->new( file => $matrix_file )
	or die " unable to open matrix file $matrix_file!\n";
unless ( $Matrix->number_columns >= 2 ) {
	die " key matrix must have at least 2 colunns!";
}
unless (@group_indexes) {
	# using remaining columns as the groups
	my $name_i = $Matrix->name_column || 0;
	my $id_i   = $Matrix->id_column || 0;
	my $n      = $Matrix->last_column;
	if ($n == 2) {
		push @group_indexes, 2;
	}
	elsif ( $name_i == 2 or $id_i == 2 ) {
		# sometimes the second column may be the name so skip it
		if ($n == 3) {
			push @group_indexes, 3;
		}
		else {
			@group_indexes = ( 3 .. $n );
		}
	}
	else {
		# use all remaining columns
		@group_indexes = ( 2 .. $n );
	}
}
printf "  Using '%s' as key and '%s' as group\n", $Matrix->name(1),
	join(', ', map { $Matrix->name($_) } @group_indexes );
$Matrix->iterate( sub {
	my $row      = shift;
	my $key      = $row->value(1);
	my $category = join '-', map { $row->value($_) } @group_indexes;
	$name2group{$key} = $category;
	$name2row{$key}   = $row->row_index;
	$group2count{$category} += 1;
} );
printf "  Loaded %d features for %d groups\n", scalar(keys %name2group),
	scalar(keys %group2count);


# Load data file
print " Loading data $data_file...\n";
my $Data = Bio::ToolBox::Data->new( file => $data_file )
	or die " unable to open data file $data_file!\n";
my $key_i = $Data->find_column( $Matrix->name(1) );
unless ($key_i) {
	die sprintf(" unable to find key index column '%s' in data file '%s'!", 
		$Matrix->name(1), $data_file);
}
$Data->iterate( sub {
	my $row = shift;
	my $name = $row->value($key_i);
	unless (exists $name2group{$name}) {
		print " Feature '$name' not in matrix list!\n";
		return;
	}

	# get category
	my $category = $name2group{$name};
	unless ( exists $group2data{$category} ) {
		$group2data{$category} = $Data->duplicate;
	}
	$group2data{$category}->add_row($row);
} );

# Identify data columns
unless (@val_indexes) {
	for my $i ( 1 .. $Data->last_column ) {
		if ( $Data->metadata($i, 'dataset') ) {
			push @val_indexes, $i;
		}
	}
	unless (@val_indexes) {
		print " No suitable data columns found! Please specify with --value option\n";
		exit 1;
	}
}

# Normalize and sort data table for each category
foreach my $category ( keys %group2data ) {

	my $GData = $group2data{$category};

	# first normalize values if requested
	if ($normlist) {
		printf " Normalizing and sorting %d features for %s...\n", $GData->number_rows,
			$category;
		$GData->iterate( sub {
			my $row = shift;
			my @norm_val = map { $row->value($_) } @norm_indexes;
			my $norm = mean(@norm_val);
			if ($norm != 0) {
				foreach my $i (@val_indexes) {
					my $v = $row->value($i);
					next unless looks_like_number($v);
					$v -= $norm;
					if ($format) {
						$v = sprintf $format, $v;
					}
					$row->value($i, $v)
				}
			}
		} );
	}
	else {
		printf " Sorting %d features for %s...\n", $GData->number_rows,
			$category;
	}
	
	# generate mean column
	my $mean_i = $GData->add_column('Mean');
	$GData->iterate( sub {
		my $row = shift;
		my @values = map { $row->value($_) } @val_indexes;
		$row->value( $mean_i, mean(@values) );
	} );
	
	# sort by descending row mean
	$GData->sort_data($mean_i, 'd');
	$GData->delete_column($mean_i);
}

# Sort the categories
my @order;
if ($sort eq 'count') {
	@order = map { $_->[1] }
				sort { $b->[0] <=> $a->[0] or $a->[1] cmp $b->[1] }
				map { [ $group2data{$_}->number_rows, $_ ] } keys %group2data;
}
elsif ($sort eq 'name') {
	@order = sort {$a cmp $b} keys %group2data;
}

# Concatenate the sorted category data into output data table
my $OutData = $Data->duplicate;
$OutData->add_comment("Sorted based on matrix key $matrix_file");
if ($normlist) {
	$OutData->add_comment( "Normalized values to reference columns $normlist" );
}
foreach my $group (@order) {
	$group2data{$group}->iterate( sub {
		my $row = shift;
		$OutData->add_row($row);
	} );
}

# Write the merged final output data
my $out_data_file = sprintf "%s%s.sorted%s", $Data->path, $Data->basename,
		$Data->extension;
my $success = $OutData->write_file($out_data_file);
if ($success) {
	print " Wrote new data file $success\n";
}

# Write individual summary files for each category
foreach my $category ( keys %group2data ) {
	my $sum_file = $group2data{$category}->summary_file(
		filename => sprintf("%s%s.sorted_%s", $Data->path, $Data->basename, $category),
		method   => 'trimmean'
	);
	if ($sum_file) {
		print " Wrote summary file $sum_file\n";
	}
}

# Write row key file
my $OutMatrix = Bio::ToolBox::Data->new( columns => [ $Matrix->list_columns ] );
# we need to translate the Name column from the matrix file to the ID column of the 
# data file to be used as row identifiers in the plot_peak_figures.R plotting script
# this is only pertinent in context of the peak calling pipeline, not necessarily 
# other situations
# use crude method to determine if a switch needs to be made
my $switch;
if ( $Matrix->name(1) eq 'Name' and $Data->name(1) eq 'Primary_ID' ) {
	$switch = 1;
	$OutMatrix->name(1, 'Primary_ID');
}
$OutData->iterate( sub {
	my $row = shift;
	my $n   = $row->value($key_i);
	my $i   = $name2row{$n};
	if ($switch) {
		my @val = $Matrix->get_row($i)->row_values();
		shift @val; # remove existing name
		$OutMatrix->add_row( [ $row->id, @val ] ); 
	}
	else {
		$OutMatrix->add_row( $Matrix->get_row($i) );
	}
} );
my $out_rows_file = sprintf "%s%s.sorted.row_groups.txt", $Data->path, 
		$Data->basename;
$success = $OutMatrix->write_file($out_rows_file);
if ($success) {
	print " Wrote new row matrix annotation file $success\n";
}


sub mean {
	my @v;
	foreach (@_) {
		push @v, $_ if looks_like_number($_);
	}
	return 0 unless @v;
	return sum0(@v)/scalar(@v);
}
