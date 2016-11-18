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
# https://github.com/tjparnell/HCI-Scripts/

use strict;
use Getopt::Long;
use List::Util qw(sum max);
use Bio::ToolBox::db_helper::bam; 

unless (@ARGV) {
	print <<END;

A script to remove excessive duplicate alignments, leaving up to X number 
of alignments at any given position, where X can be assigned by the user. 
This is in contrast to traditional duplicate removers that retain only 
one alignment per position.

The limit can be either an integer, where up to X alignments are tolerated,
or it can be a fraction (0.xx), whereupon the maximum number of alignments 
that achieve this target fraction is calculated.

Single-end aligment duplicates are checked for start position, strand, and 
calculated alignment end position to check for duplicates. Because of this, 
the numbers are slightly different than calculated by traditional duplicate 
removers.

Paired-end alignments are treated as fragments. Only properly paired 
alignments are considered; everything else is silently skipped. Fragments 
are checked for start position and fragment length (paired insertion size) 
for duplicates.

USAGE:  bam_partial_dedup.pl --max X --in in.bam --out out.bam
        bam_partial_dedup.pl --frac 0.xx --in in.bam
        bam_partial_dedup.pl -m X -i in.bam -o out.bam -d dup.bam
       
OPTIONS:
        --in    The input bam file, should be sorted and indexed
        --out   The output bam file, optional if you're just checking 
                the duplication rate for a provided fraction
        --dup   The duplicates bam file, optional if you want to keep 
                the duplicate alignments
        --max   Integer representing the maximum number of duplicates 
                at each position
        --frac  Decimal fraction representing the target duplication 
                rate in the final file. The actual maximum number will 
                be calculated prior before writing. Mutually exclusive 
                with the --max option
        --pe    Bam files contain paired-end alignments and only 
                properly paired duplicate fragments will be checked for 
                duplication. 

END
	exit;
}

### Get Options
my ($fraction, $max, $infile, $outfile, $dupfile, $paired);
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'dup=s'      => \$dupfile, # the name of the duplicates file
	'frac=f'     => \$fraction, # target fraction of duplicates
	'max=i'      => \$max, # the maximum number of alignments per position
	'pe!'        => \$paired, # treat as paired-end alignments
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

if ($max and $fraction) {
	die "max and frac options are mutually exclusive. Pick one.\n";
}
if ($max and $max !~ /^\d+$/) {
	die "unrecognized max limit '$max'. Should be integer (X)\n";
}
if ($fraction and $fraction !~ /^0\.\d+$/) {
	die "unrecognized fraction '$fraction'. Should be decimal (0.x)\n";
}
unless ($infile) {
	die "no input file provided!\n";
}




### Open bam files
# input bam file
my $sam = open_bam_db($infile) or # automatically takes care of indexing
	die "unable to read bam $infile $!";

# write header
my $header = $sam->bam->header;



### Global varaibles
my $totalCount = 0;
my %depth2count; # for generating a histogram of duplicate numbers
               # key=depth, value=number of bases with this depth
my $items = $paired ? 'properly paired fragments' : 'alignments';
	# explain what we are counting in the output




### count the reads and number at each position
if ($fraction) {
	# count alignments on each chromosome
	print " Counting $items....\n";
	for my $tid (0 .. $sam->n_targets - 1) {
		my $seq_length = $sam->target_len($tid);
	
		# prepare callback data structure
		my $data = {
			position   => -1, # a non-coordinate
			reads      => [], # temp buffer for collecting reads
		};
	
		# walk through the reads on the chromosome
		my $callback = $paired ? \&count_pe_callback : \&count_se_callback;
		$sam->bam_index->fetch($sam->bam, $tid, 0, $seq_length, $callback, $data);
	
		# check to make sure we don't leave something behind
		if (defined $data->{reads}->[0]) {
			$paired ? count_up_pe_alignments($data) : count_up_se_alignments($data);
		}
	}

	# calculate fractions
	my @depths = sort {$a <=> $b} keys %depth2count;
#	printf " Counts at depth\n%s\n", join("\n", map { sprintf " depth %s at %s positions", $_, $depth2count{$_} } sort {$a <=> $b} keys %depth2count);
	my $nondupCount = sum(values %depth2count);
	my $dupCount = $totalCount - $nondupCount;
	my $original_duplicate_rate = $dupCount / $totalCount;
	my $rate = $original_duplicate_rate;
	foreach (my $i = 1; $i <= $#depths; $i++) {
		my $depth = $depths[$i];
		my $sumCount = 0;
		foreach my $d (@depths) {
			if ($d <= $depth) {
				$sumCount += ($depth2count{$d} * $d);
			}
			else {
				$sumCount += ($depth2count{$d} * $depth);
			}
		}
		$rate = ($totalCount - $sumCount) / $totalCount;
# 		printf "  testing %s, effective duplicate rate is $rate\n", $depths[$i];
		if ($rate < $fraction) {
			$max = $depths[$i]; 
			last;
		}
	}
	printf "
 %12s total mapped $items
 %12s predicted duplicate $items
 %12s non-duplicate $items
 duplicate rate at max 1 is %.4f
 maximum duplicate depth is %d, mean depth is %.4f
 setting maximum duplicate count at $max with estimated rate of %.4f\n", 
	$totalCount, $dupCount, $nondupCount, $dupCount/$totalCount, 
	max(keys %depth2count), $totalCount / sum(values %depth2count), $rate;
}

if ($max == 1) {
	print "can't remove any more than 1!\n";
	exit;
}





### Write out new bam file with specified number of targets
unless ($outfile) {
	print " No outfile defined. Exiting.\n";
	exit;
}
print " Removing duplicates and writing new bam file....\n";

# open bam new bam file
my $outbam = write_new_bam_file($outfile) or 
	die "unable to open output bam file $outfile! $!";
	# this uses low level Bio::DB::Bam object
$outbam->header_write($header);

# open duplicate bam file if requested
my $dupbam;
if ($dupfile) {
	$dupbam = write_new_bam_file($dupfile) or 
		die "unable to open duplicate bam file $dupfile! $!";
		# this uses low level Bio::DB::Bam object
	$dupbam->header_write($header);
}

$totalCount    = 0; # reset to zero from before
my $keepCount  = 0;
my $tossCount  = 0;
for my $tid (0 .. $sam->n_targets - 1) {
	my $seq_length = $sam->target_len($tid);
	
	# prepare callback data structure
	my $data = {
		position   => -1,
		reads      => [],
		keepers    => {}, # hash of names of rev pe reads to keep
		dupkeepers => {}, # hash of names of rev pe dup reads to keep
	};
	
	# walk through the reads on the chromosome
	my $callback = $paired ? \&write_pe_callback : \&write_se_callback;
	$sam->bam_index->fetch($sam->bam, $tid, 0, $seq_length, $callback, $data);
	
	# check to make sure we don't leave something behind
	# possible with single-end, should be unlikely with paired-end
	if (defined $data->{reads}->[0]) {
		$paired ? write_out_pe_alignments($data) : write_out_se_alignments($data);
	}
	if ($paired and scalar keys %{$data->{keepers}}) {
		printf " missing reverse alignments: %s\n", 
			join(", ", keys %{$data->{keepers}});
	}
}




### Print results
printf "
 %12s total mapped $items
 %12s duplicate mapped $items discarded
 %12s mapped $items retained
 new duplicate rate %.4f
", $totalCount, $tossCount, $keepCount, $tossCount/$totalCount;




### Finish up
undef $outbam;
check_bam_index($outfile);
if ($dupfile) {
	undef $dupbam;
	check_bam_index($dupfile);
}






### single-end alignment callback for reading
sub count_se_callback {
	my ($a, $data) = @_;
	$totalCount++;
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		if ($data->{position} != -1) {
			# a real position
			count_up_se_alignments($data);
		}
		$data->{position} = $pos;
		$data->{reads} = []; # clear out existing alignments
		push @{ $data->{reads} }, $a;
	}
}

### paired-end alignment callback for reading
sub count_pe_callback {
	my ($a, $data) = @_;
	return unless $a->proper_pair; # consider only proper pair alignments
	return if $a->reversed; # count only forward alignments
	$totalCount++;
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		if ($data->{position} != -1) {
			# a real position
			count_up_pe_alignments($data);
		}
		$data->{position} = $pos;
		$data->{reads} = []; # clear out existing alignments
		push @{ $data->{reads} }, $a;
	}
}


### Count collected single-end alignments
sub count_up_se_alignments {
	my $data = shift;
	
	# split up based on end point and strand
	my %fends; # forward ends
	my %rends; # reverse ends
	foreach my $a (@{$data->{reads}}) {
		my $end = $a->calend; 
		if ($a->reversed) {
			$rends{$end} ||= [];
			push @{ $rends{$end} }, $a;
		}
		else {
			$fends{$end} ||= [];
			push @{ $fends{$end} }, $a;
		}
	}
	
	# assign count for each depth
	foreach my $pos (keys %fends) {
		$depth2count{scalar @{ $fends{$pos} }} += 1;
	}
	foreach my $pos (keys %rends) {
		$depth2count{scalar @{ $rends{$pos} }} += 1;
	}
}


### Count collected paired-end alignments
sub count_up_pe_alignments {
	my $data = shift;
	
	# split up based on reported insertion size
	my %sizes; 
	foreach my $a (@{$data->{reads}}) {
		my $s = $a->isize; 
		$sizes{$s} ||= [];
		push @{$sizes{$s}}, $a;
	}
	
	# assign count for each depth
	foreach my $s (keys %sizes) {
		$depth2count{scalar @{$sizes{$s}}} += 1;
	}
}


### single-end alignment callback for writing
sub write_se_callback {
	my ($a, $data) = @_;
	$totalCount++;
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		# write out the max number of reads
		write_out_se_alignments($data);
		
		# reset
		$data->{reads} = [];
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}


### paired-end alignment callback for writing
sub write_pe_callback {
	my ($a, $data) = @_;
	return unless $a->proper_pair; # consider only proper pair alignments
	
	# process forward reads
	$totalCount++ unless $a->reversed; # only count forward alignments
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		# write out the max number of reads
		write_out_pe_alignments($data) if defined $data->{reads}->[0];
		
		# reset
		$data->{reads} = [];
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}


### Identify and write out single-end alignments 
sub write_out_se_alignments {
	my $data = shift;
	
	# split up based on end point and strand
	my %fends; # forward ends
	my %rends; # reverse ends
	foreach my $a (@{$data->{reads}}) {
		my $end = $a->calend; 
		if ($a->reversed) {
			$rends{$end} ||= [];
			push @{ $rends{$end} }, $a;
		}
		else {
			$fends{$end} ||= [];
			push @{ $fends{$end} }, $a;
		}
	}
	
	# write forward reads
	foreach my $pos (keys %fends) {
		for (my $i = 0; $i < $max; $i++) {
			my $a = shift @{ $fends{$pos} };
			last unless $a;
			$outbam->write1($a);
			$keepCount++;
		}
		$tossCount += scalar @{ $fends{$pos} };
		if ($dupbam and scalar @{ $fends{$pos} }) {
			while (my $a = shift @{ $fends{$pos} }) {
				$dupbam->write1($a);
			}
		}
	}
	
	# write reverse reads
	foreach my $pos (keys %rends) {
		for (my $i = 0; $i < $max; $i++) {
			my $a = shift @{ $rends{$pos} };
			last unless $a;
			$outbam->write1($a);
			$keepCount++;
		}
		$tossCount += scalar @{ $rends{$pos} };
		if ($dupbam and scalar @{ $rends{$pos} }) {
			while (my $a = shift @{ $rends{$pos} }) {
				$dupbam->write1($a);
			}
		}
	}
}


### Identify and write out paired-end alignments 
sub write_out_pe_alignments {
	my $data = shift;
	
	# split up based on reported insertion size and strand
	my %f_sizes;
	my %r_sizes;
	while (my $a = shift @{$data->{reads}}) {
		my $s = $a->isize; 
		if ($a->reversed) {
			$r_sizes{$s} ||= [];
			push @{$r_sizes{$s}}, $a;
		}
		else {
			$f_sizes{$s} ||= [];
			push @{$f_sizes{$s}}, $a;
		}
	}
	
	# write out forward alignments
	foreach my $s (sort {$a <=> $b} keys %f_sizes) {
		for (my $i = 0; $i < $max; $i++) {
			my $a = shift @{ $f_sizes{$s} };
			last unless $a;
			$outbam->write1($a);
			$data->{keepers}{$a->qname} = 1; # remember name for reverse read
			$keepCount++;
		}
		$tossCount += scalar @{ $f_sizes{$s} };
		if ($dupbam and scalar @{ $f_sizes{$s} }) {
			while (my $a = shift @{ $f_sizes{$s} }) {
				$dupbam->write1($a);
				$data->{dupkeepers}{$a->qname} = 1; # remember name for rev
			}
		}
	}
	
	# write out reverse alignments
	foreach my $s (sort {$a <=> $b} keys %r_sizes) {
		while (my $a = shift @{ $r_sizes{$s} }) {
			my $name = $a->qname;
			if (exists $data->{keepers}{$name}) {
				# we have processed the forward alignment as a keeper
				# immediately write the reverse alignment
				$outbam->write1($a);
				delete $data->{keepers}{$name};
			}
			elsif (exists $data->{dupkeepers}{$name}) {
				# we have processed the forward alignment as a duplicate keeper
				# immediately write the reverse alignment
				$dupbam->write1($a);
				delete $data->{dupkeepers}{$name};
			}
		}
	}
}






sub OLD_write_pe_callback {
	my ($a, $data) = @_;
	return unless $a->proper_pair; # consider only proper pair alignments
	
	# process reverse reads
	if ($a->reversed) {
		# make sure we have written out all the forward reads before continuing
		if ($a->pos > $data->{position} and defined $data->{reads}->[0]) {
			# we've moved past the existing forward position, so write these 
			# alignments first and identify potential reverse candidates to keep
			write_out_pe_alignments($data);
			$data->{reads} = [];
		} 
		
		# we should already have processed the forward alignment 
		my $name = $a->qname;
		if (exists $data->{keepers}{$name}) {
			# we have processed the forward alignment as a keeper
			# immediately write the reverse alignment
			$outbam->write1($a);
			delete $data->{keepers}{$name};
		}
		elsif (exists $data->{dupkeepers}{$name}) {
			# we have processed the forward alignment as a duplicate keeper
			# immediately write the reverse alignment
			$dupbam->write1($a);
			delete $data->{dupkeepers}{$name};
		}
		return;
	}
	
	# process forward reads
	$totalCount++; # only count forward alignments
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		# write out the max number of reads
		write_out_pe_alignments($data) if defined $data->{reads}->[0];
		
		# reset
		$data->{reads} = [];
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}


sub OLD_write_out_pe_alignments {
	my $data = shift;
	
	# split up based on reported insertion size
	my %sizes; 
	foreach my $a (@{$data->{reads}}) {
		my $s = $a->isize; 
		$sizes{$s} ||= [];
		push @{$sizes{$s}}, $a;
	}
	
	# write out alignments
	foreach my $s (sort {$a <=> $b} keys %sizes) {
		for (my $i = 0; $i < $max; $i++) {
			my $a = shift @{ $sizes{$s} };
			last unless $a;
			$outbam->write1($a);
			$data->{keepers}{$a->qname} = 1; # remember name for reverse read
			$keepCount++;
		}
		$tossCount += scalar @{ $sizes{$s} };
		if ($dupbam and scalar @{ $sizes{$s} }) {
			while (my $a = shift @{ $sizes{$s} }) {
				$dupbam->write1($a);
				$data->{dupkeepers}{$a->qname} = 1; # remember name for rev
			}
		}
	}
}




