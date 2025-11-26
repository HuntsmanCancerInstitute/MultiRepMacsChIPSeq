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

our $VERSION = 0.4;

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

A script to collate the annotation output from `bedtools window` using
peaks and gene annotation TSS. This is intended to be used with a
window parameter to identify all neighborhood genes.

This will report any gene that overlaps (within $overlap_distance bp),
single genes to the immediate left and right within the neighborhood distance,
and all genes within the neighborhood.

Note that there will be one-to-one, one-to-many, many-to-one, and many-to-many
relationships between peaks and genes. Expect duplications.



Writes several output files:
	outfile.txt                Peak file with overlapping, left, right, and
	                             neighborhood gene annotation
	outfile.overlapping.txt    Text file of all overlapping genes with peak names
	outfile.closest.txt        Text file of the closest genes with peak names
	outfile.adjacent.txt       Text file of all overlapping, immediate left, and 
	                             immediate right genes with peak names
    outfile.neighborhood.txt   Text file of all neighborhood TSS features with
                                 gene names and peak names
    outfile.profile.txt        Text file of of peak coverage around TSS

VERSION: $VERSION

USAGE: script.pl  ...

OPTIONS:
  -i --in <file>           Provide the input file
  -o --out <file>          Provide the output file
  -a --alt                 Alternate lookup for transcript to name
  -d --distance <int>      Neighborhood distance for reporting in Kb (10)
  --overlap <int>          Gap distance to consider overlapping in bp ($overlap_distance)
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
my $OutData = Bio::ToolBox->new_data(
	qw(Coordinate Name Closest_tssID Closest_GeneID Closest_GeneName
		Overlapping_tssID Overlapping_GeneID Overlapping_GeneName
		Left_tssID LeftDistance Left_GeneID Left_GeneName
		Right_tssID RightDistance Right_GeneID Right_GeneName
		Neighborhood_GeneID Neighborhood_GeneName) 
);
my $OverlapData      = Bio::ToolBox->new_data( qw(GeneID GeneName Distance PeakID PeakName) );
my $ClosestData      = Bio::ToolBox->new_data( qw(GeneID GeneName Distance PeakID PeakName) );
my $AdjacentData     = Bio::ToolBox->new_data( qw(GeneID GeneName Distance PeakID PeakName) );
my $NeighborhoodData = Bio::ToolBox->new_data( qw(tssID GeneID GeneName Distance 
						PeakID PeakName) );
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

	# reference feature
	my $reference = Bio::ToolBox::SeqFeature->new(
		-seq_id => $items->[0]->seq_id,
		-start  => $items->[0]->start,
		-end    => $items->[0]->stop,
		-name   => $items->[0]->name,
		-id     => $items->[0]->coordinate
	);
	
	# sorted seqfeatures
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
		$leftDistance = $left->end - $reference->start;
		if ( abs($leftDistance) < $overlap_distance ) {
			# consider this overlapping
			push @overlap, $left;
			( $overlapping, $overlapping_gid, $overlapping_gnm, $overlapping_start ) =
				process_list(\@overlap);
			undef $leftDistance;
			$left = pop @left_items || undef;
		}
		if ($left) {
			if ($left->strand < 0) {
				# upstream
				$leftDistance = $left->end - $reference->start;
			}
			else {
				# downstream
				$leftDistance = $reference->start - $left->start;
			}
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
			undef $rightDistance;
		}
		if ($right) {
			if ($right->strand < 0) {
				# downstream
				$rightDistance = $right->end - $reference->end;
			}
			else {
				# upstream
				$rightDistance = $reference->end - $right->start;
			}
			( $rightName, $right_gid, $right_gnm, $right_start ) =
				process_list( [ $right ] );
		}
	}

	# closest
	my ($closest, $closest_gid, $closest_gnm, $closest_dist);
	if ($overlapping) {
		$closest      = $overlapping;
		$closest_gid  = $overlapping_gid;
		$closest_gnm  = $overlapping_gnm;
		$closest_dist = 0;
	}
	elsif ($leftName and not $rightName) {
		$closest      = $leftName;
		$closest_gid  = $left_gid;
		$closest_gnm  =  $left_gnm;
		$closest_dist = $leftDistance;
	}
	elsif ($rightName and not $leftName) {
		$closest      = $rightName;
		$closest_gid  = $right_gid;
		$closest_gnm  = $right_gnm;
		$closest_dist = $rightDistance;
	}
	elsif ($rightName and $leftName) {
		if ( abs($leftDistance) < abs($rightDistance) ) {
			$closest      = $leftName;
			$closest_gid  = $left_gid;
			$closest_gnm  = $left_gnm;
			$closest_dist = $leftDistance;
		}
		else {
			$closest      = $rightName;
			$closest_gid  = $right_gid;
			$closest_gnm  = $right_gnm;
			$closest_dist = $rightDistance;
		}
	}
	if ($closest_gid) {
		if ($closest_gid =~ /,/) {
			# sigh, need to de-list this again
			# this only happens with overlapping TSSs, never with left or right
			my @gids = split /,/, $closest_gid;
			my @gnms = split /,/, $closest_gnm;
			for my $i (0..$#gids) {
				$ClosestData->add_row( [
					$gids[$i] || '.',
					$gnms[$i] || '.',
					0,
					$reference->id,
					$reference->name
				] );
			}
		}
		else {
			$ClosestData->add_row( [
				$closest_gid || '.',
				$closest_gnm || '.',
				$closest_dist,
				$reference->id,
				$reference->name
			] );
		}
	}

	# overlapping gene output
	if (@overlap) {
		my @overlap_gid = uniqstr map { $id2alt{ $_->name }->[0] } @overlap;
		my @overlap_gnm = uniqstr map { $id2alt{ $_->name }->[1] } @overlap;
		for my $i (0 .. $#overlap_gid) {
			$OverlapData->add_row( [
				$overlap_gid[$i] || '.',
				$overlap_gnm[$i] || '.',
				0,
				$reference->id,
				$reference->name
			] );
		}	

		# also add to adjacent output
		for my $i (0 .. $#overlap_gid) {
			$AdjacentData->add_row( [
				$overlap_gid[$i] || '.',
				$overlap_gnm[$i] || '.',
				0,
				$reference->id,
				$reference->name
			] );
		}	
	}
	
	# adjacent gene output
	if ($left_gid) {
		$AdjacentData->add_row( [
			$left_gid || '.',
			$left_gnm || '.',
			$leftDistance,
			$reference->id,
			$reference->name
		] );
	}
	if ($right_gid) {
		$AdjacentData->add_row( [
			$right_gid || '.',
			$right_gnm || '.',
			$rightDistance,
			$reference->id,
			$reference->name
		] );
	}
	
	# neighorhood output - this is essentially everything
	foreach my $f (@features) {
		my $dist;
		# tss left of the peak
		if ($f->end < $reference->start) {
			if ($f->strand < 0) {
				# upstream
				$dist = $f->end - $reference->start;
			}
			else {
				# downstream
				$dist = $reference->start - $f->start;
			}
		}
		# tss right of the peak
		elsif ($f->start > $reference->end) {
			if ($f->strand < 0) {
				# downstream
				$dist = $f->end - $reference->end;
			}
			else {
				# upstream
				$dist = $reference->end - $f->start;
			}
		}
		# overlapping
		else {
			$dist = 0;
		}
		$NeighborhoodData->add_row( [
			$f->name,
			$id2alt{ $f->name }->[0] || '.',
			$id2alt{ $f->name }->[1] || '.',
			$dist,
			$reference->id,
			$reference->name
		] );
	}
	
	# peak output
	my ( undef, $neighbor_gid, $neighbor_gnm, undef ) = process_list(\@features);
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
		$neighbor_gid    || q(.),
		$neighbor_gnm    || q(.)
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

	# overlapping genes
	$OverlapData->sort_data(2, 'i');
	printf " > identified %d unique overlapping genes\n",
		scalar( uniqstr( $OverlapData->column_values(1) ) ) - 1;
	$outfile =~ s/\.txt (?:\.gz)? $//x;
	$outfile .= '.overlapping.txt';
	$w = $OverlapData->write_file($outfile);
	unless ($w) {
		print " Failed to write overlapping gene list file!\n";
	}
	undef $w;

	# closest genes
	$ClosestData->sort_data(2, 'i');
	printf " > identified %d unique closest genes\n",
		scalar( uniqstr( $ClosestData->column_values(1) ) ) - 1;
	$outfile =~ s/overlapping/closest/x;
	$w = $ClosestData->write_file($outfile);
	unless ($w) {
		print " Failed to write closest gene list file!\n";
	}
	
	# adjacent genes
	$AdjacentData->sort_data(2, 'i');
	printf " > identified %d unique adjacent genes\n",
		scalar( uniqstr( $AdjacentData->column_values(1) ) ) - 1;
	$outfile =~ s/closest/adjacent/x;
	$w = $AdjacentData->write_file($outfile);
	unless ($w) {
		print " Failed to write adjacent gene list file!\n";
	}
	
	# neighborhood genes
	$NeighborhoodData->sort_data(3, 'i');
	printf " > identified %d unique neighborhood genes\n",
		scalar( uniqstr( $NeighborhoodData->column_values(2) ) ) - 1;
	$outfile =~ s/adjacent/neighborhood/x;
	$w = $NeighborhoodData->write_file($outfile);
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

