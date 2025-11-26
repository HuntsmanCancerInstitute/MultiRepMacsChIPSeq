#!/usr/bin/env perl

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
use File::Spec;
use File::Which;
use Bio::ToolBox 1.70;
use Bio::ToolBox::utility qw(format_with_commas);
use Bio::ToolBox::SeqFeature;
use List::Util qw(uniqstr uniqint);
use Set::IntSpan::Fast;

our $VERSION = 0.2;

# user variables
my $tool = which('bedtools');
my $infile;
my $outfile;
my $alt_file;
my $distance         = 10;   # kb
my $overlap_distance = 250;  # bp
my $profile_radius;   # kb
my $help;

# global variables

### Documentation
my $docs = <<DOC;

A script to identify the overlapping genes

VERSION: $VERSION

USAGE: script.pl  ...

OPTIONS:
  -i --in <file>           Provide the input file
  -o --out <file>          Provide the output file
  -d --distance <int>      Neighborhood distance for reporting in Kb (10)
  -a --alt                 Alternate lookup for transcript to name
  -h --help                Print documentation
DOC

# Options
unless (@ARGV) {
	print $docs;
	exit;
}
GetOptions(
	'i|in=s'        => \$infile,
	'o|out=s'       => \$outfile,
	'd|distance=i'  => \$distance,
	'a|alt=s'       => \$alt_file,
	'h|help!'       => \$help,
) or die "unrecognized option!\n";

$distance *= 1000;
if ($profile_radius) {
	$profile_radius *= 1000;
}
else {
	$profile_radius = $distance;
}
unless ($infile) {
	die "no input file!";
}
unless ($outfile) {
	die "no output file!";
}
unless ($alt_file) {
	die "no alternate names file!";
}

# prepare output
my %id2alt;
my %neighborhood;
my %all_overlapping;
my %all_closest;
my $OutData = Bio::ToolBox->new_data(
	qw(Coordinate Name Closest_ID Closest_GeneID Closest_GeneName
		Overlapping_ID Overlapping_GeneID Overlapping_GeneName
		Left_ID LeftDistance Left_GeneID Left_GeneName
		Right_ID RightDistance Right_GeneID Right_GeneName
		Neighborhood_GeneID Neighborhood_GeneName) 
);
my %profile = map { $_ => 0 }
	( ( -1 * int( $profile_radius / 100 ) ) .. int( $profile_radius / 100 ) );


# process input
load_alt_names();
parse_file();
write_output_files();





############ Subroutines 

sub load_alt_names {
	my $Data = Bio::ToolBox->load_file ( file => $alt_file );
	my $trxid = $Data->find_column('transcript.?id');
	my $genid = $Data->find_column('gene.?id');
	my $gname = $Data->find_column('gene.?name');
	$Data->iterate( sub {
		my $row = shift;
		$id2alt{ $row->value($trxid) } = [ $row->value($genid), $row->value($gname) ];
	} );
}

sub parse_file {

	# Load input file
	my $Data = Bio::ToolBox->load_file( file => $infile, noheader => 1 )
		or die "unable to load file '$infile'!\n";
	printf " Loaded '$infile' with %d features\n", $Data->number_rows;
	unless ($Data->number_columns == 11) {
		die " file does not have 11 columns!";
	}

	# rename columns to make things easier
	$Data->name(1, 'Chromosome');
	$Data->name(2, 'Start0');
	$Data->name(3, 'End');
	$Data->name(4, 'Name');

	# iterate through features
	my $current = $Data->value(1,4);
	my @items;
	$Data->iterate( sub {
		my $Row = shift;
		if (@items and $current ne $Row->name) {
			process_items( \@items );
			@items = ( $Row );
			$current = $Row->name;
		}
		else {
			push @items, $Row;
		}
	} );
	
	# remaining
	if (@items) {
		process_items( \@items );
	}
}

sub process_items {
	my $items = shift;
	# items is an array reference of Bio::ToolBox::Data::Features

	# build a reference feature
	# all should have the same reference coordinates, so take the first as an example
	my $reference = Bio::ToolBox::SeqFeature->new(
		-seq_id => $items->[0]->seq_id,
		-start  => $items->[0]->start,
		-end    => $items->[0]->stop,
		-name   => $items->[0]->name,
		-id     => $items->[0]->coordinate
	);
	
	# generate sorted list of found tss seqfeatures
	my @features = 
		map  { $_->[1] }
		sort { $a->[0] cmp $b->[0] }
		map  { [$_->start, $_] }
		map  {
			Bio::ToolBox::SeqFeature->new(
				-seq_id => $_->value(6),
				-start  => $_->value(7),
				-end    => $_->value(8),
				-name   => $_->value(9),
				-strand => $_->value(11)
			)
		} @{ $items };
	
	# find overlapping 
	my @overlap = grep { $reference->overlaps($_) } @features;
	my ( $overlapping, $overlapping_gid, $overlapping_gnm, $overlapping_start ) =
		process_list(\@overlap);
	
	
	# find left item
	my @left_items = grep { $_->end < $reference->start } @features;
	my ($left, $leftName, $leftDistance, $left_gid, $left_gnm, $left_start);
	if (@left_items) {
		$left = pop @left_items; # last one should be closest
		$leftDistance = $reference->start - $left->end;
		if ( $leftDistance < $overlap_distance ) {
			# consider this overlapping
			push @overlap, $left;
			( $overlapping, $overlapping_gid, $overlapping_gnm, $overlapping_start ) =
				process_list(\@overlap);
			$leftDistance = 0;
			$left = pop @left_items || undef;
		}
		if ($left) {
			$leftDistance = $reference->start - $left->end;
			( $leftName, $left_gid, $left_gnm, $left_start ) = process_list( [ $left ] );
		}
	}

	# find right item
	my @right_items = grep { $_->start > $reference->end } @features;
	my ($rightName, $rightDistance, $right_gid, $right_gnm, $right_start);
	if (@right_items) {
		my $right = shift @right_items; # first one should be closest
		$rightDistance = $right->start - $reference->end;
		if ( $rightDistance < $overlap_distance ) {
			# consider this overlapping
			push @overlap, $right;
			( $overlapping, $overlapping_gid, $overlapping_gnm, $overlapping_start ) =
				process_list(\@overlap);
			$right = shift @right_items || undef;
			$rightDistance = 0;
		}
		if ($right) {
			$rightDistance = $right->start - $reference->end;
			( $rightName, $right_gid, $right_gnm, $right_start ) =
				process_list( [ $right ] );
		}
	}

	# closest
	my ($closest, $closest_gid, $closest_gnm);
	if ($overlapping) {
		$closest = $overlapping;
		$closest_gid = $overlapping_gid;
		$closest_gnm = $overlapping_gnm;
	}
	elsif ($leftName and not $rightName) {
		$closest = $leftName;
		$closest_gid = $left_gid;
		$closest_gnm = $left_gnm;
	}
	elsif ($rightName and not $leftName) {
		$closest = $rightName;
		$closest_gid = $right_gid;
		$closest_gnm = $right_gnm;
	}
	elsif ($rightName and $leftName) {
		if ($leftDistance < $rightDistance) {
			$closest = $leftName;
			$closest_gid = $left_gid;
			$closest_gnm = $left_gnm;
		}
		else {
			$closest = $rightName;
			$closest_gid = $right_gid;
			$closest_gnm = $right_gnm;
		}
	}
	if ($closest) {
		my @g = split /,/, $closest_gid;
		my @n = split /,/, $closest_gnm;
		for my $i (0 .. $#g) {
			if (exists $all_closest{ $g[$i] } ) {
				$all_closest{ $g[$i] }->[1]++;
			}
			else {
				$all_closest{ $g[$i] } = [ $n[$i], 1 ];
			}
		}
	}

	# neighborhood
	my @neighbors;
	my @neighbor_gid;
	my @neighbor_gnm;
	foreach my $f (@features) {
		my $left_d  = $reference->start - $f->end;
		my $right_d = $f->start - $reference->end;
		if ( $left_d > 1 and $left_d <= $distance ) {
			push @neighbors, $f;
		}
		elsif ( $right_d > 1 and $right_d <= $distance ) {
			push @neighbors, $f;
		}
		# we will also include overlaps
	}
	push @neighbors, @overlap; # overlaps count in the neighborhood too
	if (@neighbors) {
		my @neighborNms = uniqstr map { $_->name } @neighbors;
		@neighbor_gid = uniqstr map { $id2alt{$_}->[0] } @neighborNms;
		@neighbor_gnm = uniqstr map { $id2alt{$_}->[1] } @neighborNms;
		for my $i (0 .. $#neighbor_gid) {
			my $n = $neighbor_gid[$i];
			if ( exists $neighborhood{$n} ) {
				$neighborhood{$n}->[1]++;
			}
			else {
				$neighborhood{$n} = [ $neighbor_gnm[$i] , 1 ];
			}
		}
	}
	
	# overlapping genes
	if (@overlap) {
		my @overlap_gid = uniqstr map { $id2alt{ $_->name }->[0] } @overlap;
		my @overlap_gnm = uniqstr map { $id2alt{ $_->name }->[1] } @overlap;
		for my $i (0 .. $#overlap_gid) {
			my $gid = $overlap_gid[$i];
			if ( exists $all_overlapping{$gid} ) {
				$all_overlapping{$gid}->[1]++;
			}
			else {
				$all_overlapping{$gid} = [ $overlap_gnm[$i], 1];
			}
		}		
	}
	
	# output
	$OutData->add_row( [
		$reference->id,
		$reference->name,
		$closest         || q(.),
		$closest_gid     || q(.),
		$closest_gnm     || q(.),
		$overlapping     || q(.),
		$overlapping_gid || q(.),
		$overlapping_gnm || q(.),
		$leftName        || q(.),
		$leftDistance    || 0,
		$left_gid        || q(.),
		$left_gnm        || q(.),
		$rightName       || q(.),
		$rightDistance   || 0,
		$right_gid       || q(.),
		$right_gnm       || q(.),
		@neighbor_gid ? join( q(,), @neighbor_gid ) : q(.),
		@neighbor_gnm ? join( q(,), @neighbor_gnm ) : q(.)
	]);
	
	# process relative coverage
	record_coverage($reference, $overlapping_start, $left_start, $right_start);
}

sub process_list {
	my $list = shift;
	return unless @{$list};
	# name sorted
	my @sorted = sort { $a->[0] cmp $b->[0] } map { [ $_->name, $_ ] } @{$list};
	my $name  = join q(,), map { $_->[0] } @sorted;
	my $gid   = join q(,), uniqstr map { $id2alt{ $_->[0] }->[0] } @sorted;
	my $gnm   = join q(,), uniqstr map { $id2alt{ $_->[0] }->[1] } @sorted;
	my $start = join q(,), uniqint map { $_->[1]->start } @sorted;
	return ($name, $gid, $gnm, $start);
}

sub record_coverage {
	my ($reference, $overlapping_start, $left_start, $right_start) = @_;
	my $ref_span = Set::IntSpan::Fast->new;
	$ref_span->add_range( $reference->start, $reference->stop );
	if ( $overlapping_start and $overlapping_start =~ /,/ ) {
		foreach my $start ( split /,/, $overlapping_start ) {
			next unless $start;
			record_intersection($ref_span, $start);
		}
	}
	elsif ($overlapping_start) {
		record_intersection($ref_span, $overlapping_start);
	}
	if ($left_start) {
		record_intersection($ref_span, $left_start);
	}
	if ($right_start) {
		record_intersection($ref_span, $right_start);
	}
}

sub record_intersection {
	my ($ref_span, $start) = @_;
	my $s = $start - $profile_radius;
	my $e = $start + $profile_radius;
	$s = 1 if $s <= 0;
	my $set = Set::IntSpan::Fast->new;
	$set->add_range( $s, $e );
	my $overlap = $ref_span->intersection($set);
	if ($overlap) {
		# subtract the reference point to make it +/- relative
		# divide into 100 bp bins for recording
		my @bins = uniqint map { int( ( $_ - $start ) / 100 ) } $overlap->as_array;
		foreach my $i (@bins) {
			next unless exists $profile{$i};
			$profile{$i} += 1;
		}
	}
}

sub write_output_files {

	# output table
	my $w = $OutData->write_file($outfile);
	unless ($w) {
		print " Failed to write output annotation file!\n";
	}
	undef $w;

	# all overlapping gene
	my $OverlapData = Bio::ToolBox->new_data( qw(GeneID GeneName Count) );
	foreach my $gid ( sort {$a cmp $b} keys %all_overlapping ) {
		$OverlapData->add_row(
			[ $gid, $all_overlapping{$gid}->[0] || q(), $all_overlapping{$gid}->[1] ]
		);
	}
	printf " > identified %d overlapping genes\n", $OverlapData->number_rows;
	$outfile =~ s/\.txt (?:\.gz)? $//x;
	$outfile .= '.overlapping.txt';
	$w = $OverlapData->write_file($outfile);
	unless ($w) {
		print " Failed to write overlapping gene list file!\n";
	}
	undef $w;

	# all closest gene
	my $ClosestData = Bio::ToolBox->new_data( qw(GeneID GeneName Count) );
	foreach my $gid ( sort {$a cmp $b} keys %all_closest ) {
		$ClosestData->add_row(
			[ $gid, $all_closest{$gid}->[0] || q(), $all_closest{$gid}->[1] ]
		);
	}
	printf " > identified %d closest genes\n", $ClosestData->number_rows;
	$outfile =~ s/overlapping/closest/x;
	$w = $ClosestData->write_file($outfile);
	unless ($w) {
		print " Failed to write closest gene list file!\n";
	}
	
	# neighborhood genes
	my $NData = Bio::ToolBox->new_data( qw(GeneID GeneName Count) );
	foreach my $g ( sort { $a cmp $b } keys %neighborhood ) {
		$NData->add_row( [ $g, $neighborhood{$g}->[0] || q(), $neighborhood{$g}->[1] ] );
	}
	printf " > identified %d neighborhood genes\n", $NData->number_rows;
	$outfile =~ s/closest/neighborhood/x;
	$w = $NData->write_file($outfile);
	unless ($w) {
		print " Failed to write neighborhood gene list file!\n";
	}
	undef $w;

	# coverage profile
	my $outname = $outfile;
	$outname =~ s/\.neighborhood\.txt//x;
	$outname =~ s/[\._]? annot (?:ation)? //x;
	my $ProfileData = Bio::ToolBox->new_data( 'Window', 'Midpoint', $outname );
	$ProfileData->add_comment(
'Interval coverage profile relative to overlapping, leftmost, and rightmost genes'
	);
	$ProfileData->add_comment('Positions in Kb');
	foreach my $key (sort {$a <=> $b} keys %profile) {
		next if $key == 0; # weird things happen around zero
		$ProfileData->add_row( [
			sprintf("%s:%.1f", $outname, $key / 10),
			$key / 10,
			$profile{$key}
		]);
	}
	$outfile =~ s/neighborhood/profile/;
	$w = $ProfileData->write_file($outfile);
	unless ($w) {
		print " Failed to write profile summary file!\n";
	}
	undef $w;
}

