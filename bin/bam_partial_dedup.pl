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
use Bio::ToolBox::db_helper 1.50 qw(
	open_db_connection 
	low_level_bam_fetch
	$BAM_ADAPTER
);
# this can import either Bio::DB::Sam or Bio::DB::HTS depending on availability
# this script is mostly bam adapter agnostic

my $VERSION = 1.7;

unless (@ARGV) {
	print <<END;

A script to remove excessive duplicate alignments down to an acceptable 
fraction. This is in contrast to traditional duplicate removers that 
remove all duplicates, retaining only one alignment per position. 

Duplicate reads may be artificial (derived from PCR duplication due to 
very low input) or natural (multiple DNA fragments enriched at few narrow 
foci, such as a ChIP peaks). Removing all duplicates can significantly 
reduce high-enrichment peaks, but not removing duplicates can lead to 
false positives. An optimal balance is therefore desirable. This is 
most important when comparing between replicates or samples.

This script has two primary modes: 
1. Randomly subsample or remove duplicate reads to reach a target 
   duplication rate. Results in a more uniform duplicate reduction 
   across the genome, which can be more typical of true PCR duplication. 
   Set the --frac and --random options below. Can also optionally set 
   the --max option to remove extreme outliers.
2. Remove duplicate reads at positions that exceed a threshold, either 
   manually set or automatically calculated, to achieve a target 
   duplication rate. Not usually recommmended as it reduces signal at  
   top peaks without addressing low level duplication across the genome.

For ChIPSeq applications, check the duplication level of the input. For 
mammalian genomes, typically 10-20% duplication is observed in sonicated 
input. ChIP samples typically have higher duplication. Set the target 
fraction of all samples to the lowest observed duplication rate.

Single-end aligment duplicates are checked for start position, strand, and 
calculated alignment end position to check for duplicates. Because of this, 
the numbers may be slightly different than calculated by traditional duplicate 
removers.

Paired-end alignments are treated as fragments. Only properly paired 
alignments are considered; everything else is silently skipped. Fragments 
are checked for start position and fragment length (paired insertion size) 
for duplicates. Random subsampling should not result in broken pairs.

Alignment duplicate marks (bit flag 0x400) are ignored. 

USAGE:  bam_partial_dedup.pl --in in.bam
        bam_partial_dedup.pl --frac 0.xx --rand --in in.bam --out out.bam
        bam_partial_dedup.pl --m X -i in.bam -o out.bam -d dup.bam
       
OPTIONS:
    --in <file>    The input bam file, should be sorted and indexed
    --out <file>   The output bam file, optional if you're just checking 
                   the duplication rate for a provided fraction
    --dup <file>   The duplicates bam file, optional if you want to keep 
                   the duplicate alignments. Non-standard alignments 
                   are not retained.
    --random       Randomly subsamples duplicate alignments so that the  
                   final duplication rate will match target duplication 
                   rate. Must set --frac option. 
    --frac <float> Decimal fraction representing the target 
                   duplication rate in the final file. 
    --max <int>    Integer representing the maximum number of duplicates 
                   at each position
    --pe           Bam files contain paired-end alignments and only 
                   properly paired duplicate fragments will be checked for 
                   duplication. 
    --seed <int>   Provide an integer to set the random seed generator to 
                   make the subsampling consistent (non-random).

END
	exit;
}

### Get Options
my ($fraction, $max, $infile, $outfile, $dupfile, $random, $paired, $seed);
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'dup=s'      => \$dupfile, # the name of the duplicates file
	'random!'    => \$random, # flag to random downsample duplicates
	'frac=f'     => \$fraction, # target fraction of duplicates
	'max=i'      => \$max, # the maximum number of alignments per position
	'pe!'        => \$paired, # treat as paired-end alignments
	'seed=i'     => \$seed, # seed for non-random random subsampling
	'bam=s'      => \$BAM_ADAPTER, # specifically set the bam adapter, advanced!
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

if ($max and $fraction and not $random) {
	die "max and frac options are mutually exclusive without random. Pick one.\n";
}
if ($fraction and $fraction !~ /^0\.\d+$/) {
	die "unrecognized fraction '$fraction'. Should be decimal (0.x)\n";
}
unless ($infile) {
	die "no input file provided!\n";
}




### Open bam files
# input bam file
my $sam = open_db_connection($infile) or # automatically takes care of indexing
	die "unable to read bam $infile $!";

# read header and set adapter specific alignment writer
my ($header, $write_alignment);
if ($BAM_ADAPTER eq 'sam') {
	$header = $sam->bam->header;
	$write_alignment = \&write_sam_alignment;
}
elsif ($BAM_ADAPTER eq 'hts') {
	$header = $sam->hts_file->header_read;
	$write_alignment = \&write_hts_alignment;
}
else {
	die "unrecognized bam adapter $BAM_ADAPTER!";
}


### Global varaibles
# these get modified in the count_*_callbacks 
my $totalCount = 0;
my %depth2count; # for generating a histogram of duplicate numbers
               # key=depth, value=number of bases with this depth
my $items = $paired ? 'properly paired fragments' : 'alignments';
	# explain what we are counting in the output
my $write_out_alignments; # subroutine reference for writing out alignments
my $chance; # probability for subsampling

# set non-random seed
if ($seed) {
	srand($seed);
}



### count the reads and number at each position
if ($fraction > 0 or not defined $max) {
	count_alignments();
}


### Write out new bam file with specified number of targets
exit 0 unless ($outfile);
exit 0 unless ($max > 1 or $chance > 0);
print " Removing duplicates and writing new bam file....\n";

# open bam new bam file
$outfile .= '.bam' unless $outfile =~ /\.bam$/i;
my $outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
	die "unable to open output bam file $outfile! $!";
	# this uses low level Bio::DB::Bam object
	# using an unexported subroutine imported as necessary depending on bam availability
$outbam->header_write($header);

# open duplicate bam file if requested
my $dupbam;
if ($dupfile) {
	$dupfile .= '.bam' unless $dupfile =~ /\.bam$/i;
	$dupbam = Bio::ToolBox::db_helper::write_new_bam_file($dupfile) or 
		die "unable to open duplicate bam file $dupfile! $!";
		# this uses low level Bio::DB::Bam object
		# using an unexported subroutine imported as necessary depending on bam availability
	$dupbam->header_write($header);
}

# set callbacks and subroutines
my $callback = $paired ? \&write_pe_callback : \&write_se_callback;
if ($random and defined $chance) {
	$write_out_alignments = $paired ? \&write_out_random_pe_alignments : 
		\&write_out_random_se_alignments;
}
else {
	$write_out_alignments = $paired ? \&write_out_max_pe_alignments : 
		\&write_out_max_se_alignments;
}

# set counters
$totalCount        = 0; # reset to zero from before
my $positionCount  = 0;
my $duplicateCount = 0;
my $tossCount      = 0;
my $unpairedCount  = 0;

# walk through bam file again
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
	low_level_bam_fetch($sam, $tid, 0, $seq_length, $callback, $data);
	
	# check to make sure we don't leave something behind
	# possible with single-end, should be unlikely with paired-end
	if (defined $data->{reads}->[0]) {
		&$write_out_alignments($data);
	}
	if ($paired and scalar keys %{$data->{keepers}}) {
		$unpairedCount += scalar(keys %{$data->{keepers}});
	}
}




### Print results
printf "  Total mapped: %22d 
  Non-duplicate count: %15d
  Retained duplicate count: %10d 
  Removed duplicate count: %11d
  Current duplication rate: %10.4f\n",
	$totalCount, $positionCount, $duplicateCount, $tossCount, 
	$duplicateCount/($duplicateCount + $positionCount);
if ($paired and $unpairedCount) {
	printf "  Unpaired count: %20d\n", $unpairedCount;
}



### Finish up
printf " Wrote %d alignments to $outfile\n", $paired ? 
	($positionCount + $duplicateCount) * 2 : $positionCount + $duplicateCount;
printf " Wrote %d alignments to $dupfile\n", $paired ? $tossCount * 2 : $tossCount 
	if $dupfile;
exit; # bam files should automatically be closed
	# there's a bug in Bio::DB::HTS v2.9 with htslib 1.4.1 that causes duplicate 
	# alignments to be written to dupbam for some inexplicable reason when I attempt 
	# to undef the adapter and index the files.
	# I also get errors with Bio::DB::Sam, again only with dupbam
	# so for now do not index




### Count alignments and calculate fractions
sub count_alignments {
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
		low_level_bam_fetch($sam, $tid, 0, $seq_length, $callback, $data);
	
		# check to make sure we don't leave something behind
		if (defined $data->{reads}->[0]) {
			$paired ? count_up_pe_alignments($data) : count_up_se_alignments($data);
		}
	}
	
	# report duplication statistics
		# remember that depth2count is simply depth and number of positions at that depth
		# not the number of alignments
		# to get total alignments sum( map {$_ * $depth2count{$_}} keys %depth2count )
	my @depths = sort {$a <=> $b} keys %depth2count;
	my $nondupCount = sum(values %depth2count);# assumes one unique at every single position
	my $dupCount = $totalCount - $nondupCount; # essentially count of positions with 2+ alignments
	my $dupRate = $dupCount / $totalCount;
	my $maxObserved = max(keys %depth2count);
	# right justify the numbers in the printf to make it look pretty
	printf "  Total mapped: %22d
  Non-duplicate count: %15d
  Duplicate count: %19d
  Duplication rate: %18.4f
  Maximum position depth: %12d
  Mean position depth: %15.4f\n", 
		$totalCount, $nondupCount, $dupCount, $dupRate, 
		$maxObserved, $totalCount / $nondupCount;
	
	# check if we need to continue
	if ($fraction and $dupRate <= $fraction) {
		print " Actual duplication rate is less than target fraction. No de-duplication necessary.\n";
		if ($max and $max >= 1) {
			if ($maxObserved <= $max) {
				print " Maximum observed depth is less than maximum allowed. No de-duplication necessary.\n";
				exit 0;
			}
			else {
				# we will continue with max cutoff
				print " Maximum observed depth is above maximum allowed.\n";
			}
		}
		return;
	}
	if ($max and $max >= 1 and $max >= $maxObserved) {
		print " Maximum observed depth is less than maximum allowed. No maximum cutoff necessary.\n";
		undef $max;
	}
	
	# calculate fractions
	if ($fraction and not $random) {
		my $rate;
		foreach (my $i = 0; $i <= $#depths; $i++) {
			my $depth = $depths[$i];
			my $sumCount = 0; # count of acceptable reads at given depth
			foreach my $d (@depths) {
				if ($d <= $depth) {
					# we can take up to this depth
					$sumCount += ($depth2count{$d} * $d);
				}
				else {
					# throw away the excess
					$sumCount += ($depth2count{$d} * $depth);
				}
			}
			$rate = ($totalCount - $sumCount) / $totalCount;
			printf "  At maximum depth $depth: %d kept, effective duplicate rate %.4f\n", 
				$sumCount, $rate;
			if ($rate < $fraction) {
				$max = $depths[$i]; 
				last;
			}
		}
		print " Setting maximum allowed duplicates to $max\n";
	}
	elsif ($fraction and $random) {
		
		# this is the number of duplicates available to be subsampled
		my $availableDups = $dupCount;
		if ($max) {
			# adjust for the number of duplicates to be tossed for exceeding max cutoff
			foreach my $d (keys %depth2count) {
				$availableDups -= (($d - $max) * $depth2count{$d}) 
					if ($max > 1 and $d > $max);
			}
		}
		
		# this is what the new final total should be
		my $newTotal = int($nondupCount / (1 - $fraction) );
		printf "  Expected new total: %16d\n", $newTotal;
		
		# this is how many can be duplicates
		my $allowedDups = int($newTotal * $fraction);
		printf "  Expected duplicates kept: %10d\n", $allowedDups;
		
		# this is the chance for being kept
		$chance = sprintf("%.8f", $allowedDups / $availableDups); 
		printf "  Probability of keeping: %12s\n", $chance;
	}
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
	return unless $a->tid == $a->mtid;
	return if $a->reversed; # count only forward alignments
	return unless $a->mreversed;
	return if $a->isize < 0; # wierd RF pair, a F alignment should only have + isize
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
		# write out the reads
		&$write_out_alignments($data);
		
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
	return unless $a->tid == $a->mtid;
	
	# check pair orientation
	if ($a->reversed) {
		return if $a->isize > 0;
		return if $a->mreversed;
	}
	else {
		return if $a->isize < 0;
		return unless $a->mreversed;
	}
	
	# process forward reads
	$totalCount++ unless $a->reversed; # only count forward alignments
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		# write out the reads
		&$write_out_alignments($data) if defined $data->{reads}->[0];
		
		# reset
		$data->{reads} = [];
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}


### Identify and write out single-end alignments 
sub write_out_max_se_alignments {
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
		$positionCount++;
		
		# always write one alignment
		my $a1 = shift @{ $fends{$pos} };
		&$write_alignment($outbam, $a1);
		
		# write remaining alignments up to max
		for (my $i = 0; $i < $max; $i++) {
			my $a = shift @{ $fends{$pos} };
			last unless $a;
			&$write_alignment($outbam, $a);
			$duplicateCount++;
		}
		$tossCount += scalar @{ $fends{$pos} };
		if ($dupbam and scalar @{ $fends{$pos} }) {
			while (my $a = shift @{ $fends{$pos} }) {
				&$write_alignment($dupbam, $a);
			}
		}
	}
	
	# write reverse reads
	foreach my $pos (keys %rends) {
		$positionCount++;
		
		# always write one alignment
		my $a1 = shift @{ $rends{$pos} };
		&$write_alignment($outbam, $a1);
		
		# write remaining alignments up to max
		for (my $i = 0; $i < $max; $i++) {
			my $a = shift @{ $rends{$pos} };
			last unless $a;
			&$write_alignment($outbam, $a);
			$duplicateCount++;
		}
		$tossCount += scalar @{ $rends{$pos} };
		if ($dupbam and scalar @{ $rends{$pos} }) {
			while (my $a = shift @{ $rends{$pos} }) {
				&$write_alignment($dupbam, $a);
			}
		}
	}
}


### Identify and write out paired-end alignments 
sub write_out_max_pe_alignments {
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
		$positionCount++;
		
		# always write one alignment
		my $a1 = shift @{ $f_sizes{$s} };
		&$write_alignment($outbam, $a1);
		$data->{keepers}{$a1->qname} = 1; # remember name for reverse read
		
		# write remaining up to max
		for (my $i = 0; $i < $max; $i++) {
			my $a = shift @{ $f_sizes{$s} };
			last unless $a;
			&$write_alignment($outbam, $a);
			$data->{keepers}{$a->qname} = 1; # remember name for reverse read
			$duplicateCount++;
		}
		$tossCount += scalar @{ $f_sizes{$s} };
		if ($dupbam and scalar @{ $f_sizes{$s} }) {
			while (my $a = shift @{ $f_sizes{$s} }) {
				&$write_alignment($dupbam, $a);
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
				&$write_alignment($outbam, $a);
				delete $data->{keepers}{$name};
			}
			elsif (exists $data->{dupkeepers}{$name}) {
				# we have processed the forward alignment as a duplicate keeper
				# immediately write the reverse alignment
				&$write_alignment($dupbam, $a);
				delete $data->{dupkeepers}{$name};
			}
		}
	}
}


sub write_out_random_se_alignments {
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
		if (scalar @{ $fends{$pos} } == 1) {
			# there's only one alignment, automatic keep
			$positionCount++;
			&$write_alignment($outbam, $fends{$pos}->[0]);
		}
		else {
			random_se_alignment_writer($fends{$pos});
		}
	}
	
	# write reverse reads
	foreach my $pos (keys %rends) {
		if (scalar @{ $rends{$pos} } == 1) {
			# there's only one alignment, automatic keep
			$positionCount++;
			&$write_alignment($outbam, $rends{$pos}->[0]);
		}
		else {
			random_se_alignment_writer($rends{$pos});
		}
	}

}


sub random_se_alignment_writer {
	my $alignments = shift;
	
	# always write the first one
	$positionCount++;
	&$write_alignment($outbam, shift @{$alignments});
	
	# duplicates get to roll the dice
	my @keep;
	my @toss;
	while (my $a = shift @{$alignments}) {
		if (rand(1) < $chance) {
			push @keep, $a;
		}
		else {
			push @toss, $a;
		}
	}
	
	# cut off above maximum
	if ($max and scalar @keep > $max) {
		my @removed = splice(@keep, $max - 1);
		push @toss, @removed;
	}
	
	# write out keepers
	$duplicateCount += scalar(@keep);
	foreach (@keep) {
		&$write_alignment($outbam, $_);
	}
	
	# record losers
	$tossCount += scalar(@toss);
	if ($dupbam) {
		foreach (@toss) {
			&$write_alignment($dupbam, $_);
		}
	}
}

sub write_out_random_pe_alignments {
	my $data = shift;
	
	# split up based on reported insertion size and strand
	my %f_sizes;
	my %r_sizes;
	while (my $a = shift @{$data->{reads}}) {
		my $s = $a->isize; # positive for F reads, negative for R reads
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
		if (scalar @{ $f_sizes{$s} } == 1) {
			# there's only one fragment, automatic keep
			$positionCount++;
			my $a = $f_sizes{$s}->[0];
			&$write_alignment($outbam, $a);
			$data->{keepers}{$a->qname} = 1; # remember name for reverse read
		}
		else {
			# always write the first one
			$positionCount++;
			my $a = shift @{ $f_sizes{$s} };
			&$write_alignment($outbam, $a);
			$data->{keepers}{$a->qname} = 1; # remember name for reverse read
			
			# duplicates get to roll the dice
			my @keep;
			my @toss;
			while (my $a = shift @{ $f_sizes{$s} }) {
				if (rand(1) < $chance) {
					push @keep, $a;
				}
				else {
					push @toss, $a;
				}
			}
	
			# cut off above maximum
			if ($max and scalar @keep > $max) {
				my @removed = splice(@keep, $max - 1);
				push @toss, @removed;
			}
	
			# write out keepers
			$duplicateCount += scalar(@keep);
			foreach (@keep) {
				&$write_alignment($outbam, $_);
				$data->{keepers}{$_->qname} = 1; # remember name to write reverse read
			}
	
			# record losers
			$tossCount += scalar(@toss);
			if ($dupbam) {
				foreach (@toss) {
					&$write_alignment($dupbam, $_);
					$data->{dupkeepers}{$_->qname} = 1; # remember name 
				}
			}
		}
	}
	
	# write out reverse alignments
	foreach my $s (sort {$a <=> $b} keys %r_sizes) {
		# we are sorting negative sizes, so the biggest, ie those with mate furthest 
		# leftward or largest negative should go first
		while (my $a = shift @{ $r_sizes{$s} }) {
			my $name = $a->qname;
			if (exists $data->{keepers}{$name}) {
				# we have processed the forward alignment as a keeper
				# immediately write the reverse alignment
				&$write_alignment($outbam, $a);
				delete $data->{keepers}{$name};
			}
			elsif (exists $data->{dupkeepers}{$name}) {
				# we have processed the forward alignment as a duplicate keeper
				# immediately write the reverse alignment
				&$write_alignment($dupbam, $a);
				delete $data->{dupkeepers}{$name};
			}
		}
	}
}


sub write_sam_alignment {
	# wrapper for writing Bio::DB::Sam alignments
	return $_[0]->write1($_[1]);
}

sub write_hts_alignment {
	# wrapper for writing Bio::DB::HTS alignments
	return $_[0]->write1($header, $_[1]);
}



