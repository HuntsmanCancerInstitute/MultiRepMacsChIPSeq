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
use Bio::ToolBox 2.0;
use Bio::ToolBox::utility qw(format_with_commas);
use Bio::ToolBox::SeqFeature;
use List::Util qw(uniqstr uniqint any);
use Set::IntSpan::Fast;
use Data::Dumper;

our $VERSION = 0.6;

# user variables
my $tool = which('bedtools');
my $infile;
my $outfile;
my $tss_file;
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
  -i --in <file>           Provide the input peak file (bed, narrowPeak)
  -o --out <file>          Provide the output file basename
  -t --tss <file>          Provide the TSS data file
  -d --distance <int>      Neighborhood distance for reporting in Kb (10)
  --overlap <int>          Gap distance to consider overlapping in bp ($overlap_distance)
  -b --bedtools <path>     Path to bedtools ($tool)
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
	't|tss=s'       => \$tss_file,
	'b|bedtools=s'  => \$tool,
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
unless ($tss_file) {
	die "no alternate names file!";
}
unless ($tool) {
	die "no bedtools application path provided or found!";
}

# prepare output
my %chroms;
my %id2alt;
my $OutData = Bio::ToolBox->new_data( qw(Coordinate Name 
		Overlapping_GeneID Overlapping_GeneName
		Closest_GeneID Closest_GeneName Closest_Distance
		Left_GeneID Left_GeneName LeftDistance
		Right_GeneID Right_GeneName RightDistance
		Neighborhood_GeneID Neighborhood_GeneName) 
);
my %gene2peak;		# a hash of geneIDs to array of peaks
					# peaks [geneName, overlapping, closest, adjacent, neighborhood]
my @columns = qw(GeneID GeneName ClosestDistance Transcripts PeakID PeakName);
my $OverlapData      = Bio::ToolBox->new_data( @columns );
my $ClosestData      = Bio::ToolBox->new_data( @columns );
my $AdjacentData     = Bio::ToolBox->new_data( @columns );
my $NeighborhoodData = Bio::ToolBox->new_data( @columns );
my %profile = map { $_ => 0 }
	( -1 * downsize($profile_radius) .. downsize( $profile_radius ) );
my $GeneData = Bio::ToolBox->new_data( qw(GeneID GeneName Overlapping Closest
		Adjacent Neighborhood) );

# temporary output files
my $tss_bed_file  = sprintf("tss.%s.bed", $PID);
my $peak_bed_file = sprintf("peak.%s.bed", $PID);
my $result_file   = sprintf("result.%s.txt", $PID);


# main execution functions
load_peak_data();
load_tss_data();
run_intersection();
parse_intersection_file();
write_output_files();
unlink($tss_bed_file, $peak_bed_file, $result_file);




############ Subroutines 

sub load_peak_data {
	my $Data = Bio::ToolBox->load_file( file => $infile ) or
		die 'unable to load input peak file!';
	printf " Loaded peak file with %d features\n", $Data->number_rows;
	my $Peak = Bio::ToolBox->new_data( qw( Chromosome Start0 End Name ) );
	$Data->iterate( sub {
		my $row  = shift;
		$chroms{ $row->seq_id } += 1;
		my @data = ( $row->seq_id, $row->start - 1, $row->end, $row->name );
		$Peak->add_row( \@data );
	} );
	$Peak->gsort_data;
	$Peak->write_file($peak_bed_file);
}

sub load_tss_data {
	my $Data  = Bio::ToolBox->load_file ( file => $tss_file ) or
		die 'unable to load TSS annotation file!';
	printf " Loaded TSS annotation file with %d features\n", $Data->number_rows;
	my $chr   = $Data->chromo_column;
	my $start = $Data->start_column;
	my $strnd = $Data->strand_column;
	my $trxid = $Data->find_column('transcript.?id');
	my $genid = $Data->find_column('gene.?id');
	my $gname = $Data->find_column('gene.?name');
	my $TSS   = Bio::ToolBox->new_data(qw( Chromosome Start0 End Name Score Strand ) );
	$TSS->bed(6);
	$TSS->interbase(1);
	$Data->iterate( sub {
		my $row = shift;
		
		# check the chromosome first
		# only take those that we are interested in
		# also check for chr prefix
		my $chr = $row->seq_id;
		unless (exists $chroms{$chr} ) {
			my $alt;
			if ( $chr =~ /^chr (.+) $/xi ) {
				$alt = $1;
			}
			else {
				$alt = 'chr' . $chr;
			}
			if ( exists $chroms{$alt} ) {
				$chr = $alt;
			}
			else {
				next;
			}
		}
		$id2alt{ $row->value($trxid) } = [ $row->value($genid), $row->value($gname) ];
		$TSS->add_row( [
			$chr,
			$row->start - 1,
			$row->start,
			$row->value($trxid),
			1,
			$row->strand
		] );
	} );
	$TSS->gsort_data;
	$TSS->write_file($tss_bed_file);
}

sub run_intersection {
	my $command = sprintf "%s windows -a %s -b %s -w %s > %s", $tool, $peak_bed_file,
		$tss_bed_file, $distance, $result_file;
	print " Intersecting Peaks with TSS...\n";
	system($command) == 0 or
		die " execution '$command' failed!\n";
}

sub parse_intersection_file {

	# Load input file
	my $Data = Bio::ToolBox->load_file( file => $result_file, noheader => 1 )
		or die "unable to load file '$result_file'!\n";
	printf " Loaded result file with %d features\n", $Data->number_rows;

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

	# generate the peak reference feature
	# these should all be the same peak, so just take the first one
	my $reference = Bio::ToolBox::SeqFeature->new(
		-seq_id => $items->[0]->seq_id,
		-start  => $items->[0]->start,
		-end    => $items->[0]->stop,
		-name   => $items->[0]->name,
		-id     => $items->[0]->coordinate
	);
	
	# generate seqfeatures for the transcripts and sort them by position
	my @features = 
		map  { $_->[1] }
		sort { $a->[0] <=> $b->[0] }
		map  { [
			$_->value(7),
			Bio::ToolBox::SeqFeature->new(
				-seq_id     => $_->value(5),
				-start      => $_->value(6) + 1,
				-end        => $_->value(7),
				-name       => $_->value(8),
				-strand     => $_->value(9),
				-attributes => {
					gid   => $id2alt{ $_->value(8) }->[0],
					gname => $id2alt{ $_->value(8) }->[1],
				}
			)
		] }
		@{ $items };
	
	# calculate and store distance of each peak relative to each feature
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
		$f->add_tag_value('dist', $dist);
	}
	
	
	# find overlapping 
	my @overlap = grep { $reference->overlaps($_) } @features;
	
	
	# find left item
	my @left_items  = grep { $_->end < $reference->start } @features;
	my ($left, $leftDistance);
	if (@left_items) {
		$left = pop @left_items; # last one should be closest
		$leftDistance = $left->end - $reference->start;
		if ( abs($leftDistance) < $overlap_distance ) {
			# consider this overlapping
			push @overlap, $left;
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
		}
	}

	# find right item
	my @right_items = grep { $_->start > $reference->end } @features;
	my ($right, $rightDistance);
	if (@right_items) {
		$right = shift @right_items; # first one should be closest
		$rightDistance = $right->start - $reference->end;
		if ( $rightDistance < $overlap_distance ) {
			# consider this overlapping
			push @overlap, $right;
			undef $rightDistance;
			$right = shift @right_items || undef;
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
		}
	}

	# write overlap
	if (@overlap) {
		my @info = extract_info( \@overlap );
		foreach my $item (@info) {
			push @{$item}, $reference->primary_id, $reference->name;
			$OverlapData->add_row($item);
			if (exists $gene2peak{ $item->[0] } ) {
				push @{ $gene2peak{ $item->[0] }->[1] }, $reference->name;
			}
			else {
				$gene2peak{ $item->[0] } = [
					$item->[1],
					[ $reference->name ],
					[],
					[],
					[]
				];
			}
		}
	}


	# closest list includes all overlapping, and whichever left or right is closest
	my @closest;
	my $closestDistance;
	if (@overlap) {
		push @closest, @overlap;
		$closestDistance = 0;
	}
	elsif ($left and not $right) {
		push @closest, $left;
		$closestDistance = $leftDistance;
	}
	elsif ($right and not $left) {
		push @closest, $right;
		$closestDistance = $rightDistance;
	}
	elsif ($right and $left) {
		if ( abs($leftDistance) < abs($rightDistance) ) {
			push @closest, $left;
			$closestDistance = $leftDistance;
		}
		else {
			push @closest, $right;
			$closestDistance = $rightDistance;
		}
	}
	if (@closest) {
		my @info = extract_info( \@closest );
		foreach my $item (@info) {
			push @{$item}, $reference->primary_id, $reference->name;
			$ClosestData->add_row($item);
			if (exists $gene2peak{ $item->[0] } ) {
				push @{ $gene2peak{ $item->[0] }->[2] }, $reference->name;
			}
			else {
				$gene2peak{ $item->[0] } = [
					$item->[1],
					[],
					[ $reference->name ],
					[],
					[]
				];
			}
		}
	}


	# adjacent list includes all overlapping and both left and right
	my @adjacent;
	if (@overlap) {
		push @adjacent, @overlap;
	}
	if ($left) {
		push @adjacent, $left;
	}
	if ($right) {
		push @adjacent, $right;
	}
	if (@adjacent) {
		my @info = extract_info( \@adjacent );
		foreach my $item (@info) {
			push @{$item}, $reference->primary_id, $reference->name;
			$AdjacentData->add_row($item);
			if (exists $gene2peak{ $item->[0] } ) {
				push @{ $gene2peak{ $item->[0] }->[3] }, $reference->name;
			}
			else {
				$gene2peak{ $item->[0] } = [
					$item->[1],
					[],
					[],
					[ $reference->name ],
					[]
				];
			}
		}
	}


	# neighorhood output - this is essentially everything
	{
		my @info = extract_info( \@features );
		foreach my $item (@info) {
			push @{$item}, $reference->primary_id, $reference->name;
			$NeighborhoodData->add_row($item);
			if (exists $gene2peak{ $item->[0] } ) {
				push @{ $gene2peak{ $item->[0] }->[4] }, $reference->name;
			}
			else {
				$gene2peak{ $item->[0] } = [
					$item->[1],
					[],
					[],
					[],
					[ $reference->name ]
				];
			}
		}
	}


	# Peak output
	$OutData->add_row( [
		$reference->id,
		$reference->name,
		@overlap  ? join(',', map { $_->get_tag_values('gid') } @overlap ) : '.',
		@overlap  ? join(',', map { $_->get_tag_values('gname') } @overlap ) : '.',
		@closest  ? join(',', map { $_->get_tag_values('gid') } @closest ) : '.',
		@closest  ? join(',', map { $_->get_tag_values('gname') } @closest ) : '.',
		@closest  ? $closestDistance : '.',
		$left     ? $left->get_tag_values('gid') : '.',
		$left     ? $left->get_tag_values('gname') : '.',
		$left     ? $leftDistance : '.',
		$right    ? $right->get_tag_values('gid') : '.',
		$right    ? $right->get_tag_values('gname') : '.',
		$right    ? $rightDistance : '.',
		@features ? join(',', map { $_->get_tag_values('gid') } @features ) : '.',
		@features ? join(',', map { $_->get_tag_values('gname') } @features ) : '.',
	] );


	# process relative coverage
	if (@closest) {
		record_coverage( $reference, $closest[0] );
	}
}

sub extract_info {

	# this collapses the list of found features into unique genes
	# returns an array of genes consisting of (GeneID, GeneName, closest relative
	# distance to a peak, and comma-delimited list of TranscriptIDs)
	my $list = shift;

	# only one item, very easy, just return it
	if ( scalar @{$list} == 1 ) {
		my $f = $list->[0];
		return [
			$f->get_tag_values('gid'),
			$f->get_tag_values('gname'),
			$f->get_tag_values('dist'),
			$f->name
		];
	}

	# otherwise process list
	# first sort by gene id then by absolute distance to reference
	my @sorted =
		map { $_->[2] }
		sort { $a->[0] cmp $b->[0] or abs( $a->[1] ) <=> abs( $b->[1] ) }
		map { [
				$_->get_tag_values('gid'),
				$_->get_tag_values('dist'),
				$_
		] }
		@{$list};

	# then extract information, seeding with the first item
	my @info;
	push @info, [
		$sorted[0]->get_tag_values('gid'),
		$sorted[0]->get_tag_values('gname'),
		$sorted[0]->get_tag_values('dist'),
		[ $sorted[0]->name ]
	];
	for my $i (1..$#sorted) {
		my $gid = $sorted[$i]->get_tag_values('gid');
		if ($info[-1]->[0] eq $gid) {
			push @{ $info[-1][3] }, $sorted[$i]->name;
		}
		else {
			push @info, [
				$gid,
				$sorted[$i]->get_tag_values('gname'),
				$sorted[$i]->get_tag_values('dist'),
				[ $sorted[$i]->name ]
			];
		}
	}

	# finally collapse the transcript name array
	foreach my $i (@info) {
		$i->[3] = join( ',', @{ $i->[3] } );
	}
	return @info;
}


sub record_coverage {
	
	# record the coverage of the peak relative to the closest gene promoter
	my ($reference, $closest) = @_;
	
	# initialize span sets
	# we divide all coordinates into bins of 100 bp to reduce resolution
	my $ref_span = Set::IntSpan::Fast->new;
	$ref_span->add_range( downsize($reference->start), downsize($reference->stop) );
	my $tss      = downsize( $closest->start );
	my $rad      = downsize( $profile_radius );
	my $start    = $tss - $rad;
	$start       = 1 if $start <= 0;
	my $tss_span = Set::IntSpan::Fast->new;
	$tss_span->add_range( $start, $tss + $rad );

	# record intersection
	# subtract the tss point to make it relative to the gene
	my $overlap = $tss_span->intersection($ref_span);
	my @bins    = map { $_ - $tss } $overlap->as_array;
	if ( $closest->strand < 0 ) {
		# reverse if promoter is reverse strand
		@bins = map { $_ * -1 } @bins;
	}
	foreach my $b (@bins) {
		next unless exists $profile{$b};
		$profile{$b} += 1;
	}
}

sub downsize {
	my $v = shift;
	return sprintf "%.0f", ( $v / 100 );
}

sub write_output_files {

	# output peak table
	my $w = $OutData->write_file($outfile);
	unless ($w) {
		print " Failed to write output peak annotation file!\n";
	}
	undef $w;

	# output gene table
	foreach my $gid ( sort { $a cmp $b } keys %gene2peak ) {
		my $peak = $gene2peak{$gid};
		my @data = ( $gid, $peak->[0] );
		for my $i (1..4) {
			if ( scalar @{ $peak->[$i] } ) {
				push @data, join(',', @{ $peak->[$i] } );
			}
			else {
				push @data, '.';
			}
		}
		$GeneData->add_row( \@data );
	}
	$outfile =~ s/\.txt (?:\.gz)? $//x;
	$outfile .= '.genes.txt';
	$w = $GeneData->write_file($outfile);
	unless ($w) {
		print " Failed to write output gene annotation file!\n";
	}
	undef $w;
	
	# overlapping genes
	$OverlapData->delete_column(3);   # do not distance column, all zeroes
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
	$NeighborhoodData->sort_data(2, 'i');
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

=cut

I want peak-centered file:
	with overlapping, closest, left, right, and neighborhood genes
	genes are comma-delimited items

I want gene-centered file:
	with overlapping, closest, adjacent, and neighborhood peaks
	peaks are comma-delimited items




