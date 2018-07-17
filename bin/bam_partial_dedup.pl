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
use File::Which;
use List::Util qw(sum max);
use Bio::ToolBox::db_helper 1.50 qw(
	open_db_connection 
	low_level_bam_fetch
	$BAM_ADAPTER
);
# this can import either Bio::DB::Sam or Bio::DB::HTS depending on availability
# this script is mostly bam adapter agnostic
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};

my $VERSION = 1.9;

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

Since repetitive and high copy genomic regions are a big source of duplicate 
alignments, these regions can be entirely skipped by providing a file with 
recognizable coordinates. Any alignments overlapping these intervals are skipped 
in both counting and writing. 

USAGE:  bam_partial_dedup.pl --in in.bam
        bam_partial_dedup.pl --frac 0.xx --rand --in in.bam --out out.bam
        bam_partial_dedup.pl --m X -i in.bam -o out.bam
       
OPTIONS:
    --in <file>    The input bam file, should be sorted and indexed
    --out <file>   The output bam file containing unique and retained 
                   duplicates; optional if you're just checking the 
                   duplication rate.
    --mark         Write all alignments to output and mark duplicates 
                   with flag bit 0x400.
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
    --blacklist <file> Provide a bed/gff/text file of repeat regions to skip
    --chrskip <regex>  
    --seed <int>   Provide an integer to set the random seed generator to 
                   make the subsampling consistent (non-random).
    --cpu <int>    Specify the number of threads to use (4) 

END
	exit;
}

### Get Options
my ($fraction, $max, $infile, $outfile, $mark, $random, $paired, $chr_exclude, 
	$black_list, $seed, $cpu, $no_sam);
my @program_options = @ARGV;
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'mark!'      => \$mark, # mark the duplicates
	'random!'    => \$random, # flag to random downsample duplicates
	'frac=f'     => \$fraction, # target fraction of duplicates
	'max=i'      => \$max, # the maximum number of alignments per position
	'pe!'        => \$paired, # treat as paired-end alignments
	'chrskip=s'  => \$chr_exclude, # skip chromosomes
	'blacklist=s' => \$black_list, # file of high copy repeat regions to avoid
	'seed=i'     => \$seed, # seed for non-random random subsampling
	'cpu=i'      => \$cpu, # number of cpu cores to use
	'bam=s'      => \$BAM_ADAPTER, # specifically set the bam adapter, advanced!
	'nosam!'     => \$no_sam, # avoid using external sam adapter, advanced!
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
if ($parallel and not defined $cpu) {
	$cpu = 4;
}
$outfile .= '.bam' unless $outfile =~ /\.bam$/i;



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
my $items = $paired ? 'properly paired fragments' : 'alignments';
	# explain what we are counting in the output
my $write_out_alignments; # subroutine reference for writing out alignments
my $callback;
my $chance; # probability for subsampling
# my $htext = $header->text;
# $htext .= sprintf("\@PG\tID:bam_partial_dedup\tVN:%s\tCL:%s\n", 
# 	$VERSION, join(' ', $0, @program_options));
my $black_list_hash = process_black_list();

# chromosome list
my @tid_list;
for my $tid (0 .. $sam->n_targets - 1) {
	my $chr = $sam->target_name($tid);
	if ($chr_exclude and $chr =~ /$chr_exclude/i) {
		next;
	}
	push @tid_list, $tid;
}



### count the reads and number at each position
if ($fraction > 0 or not defined $max) {
	count_alignments();
}



### Write out new bam file with specified number of targets
exit 0 unless ($outfile);
exit 0 unless ($max > 0 or $chance > 0);
my ($counts, $outbam) = deduplicate();



### Print results and finish
printf "  Total mapped: %22d 
  Non-duplicate count: %15d
  Retained duplicate count: %10d 
  Removed duplicate count: %11d
  Current duplication rate: %10.4f\n",
	$counts->{total}, $counts->{nondup}, $counts->{duplicate}, $counts->{toss}, 
	$counts->{duplicate}/($counts->{duplicate} + $counts->{nondup});
if ($paired and $counts->{unpaired}) {
	printf "  Unpaired count: %20d\n", $counts->{unpaired};
}
printf " Wrote %d alignments to $outfile\n", $paired ? 
	($counts->{nondup} + $counts->{duplicate}) * 2 : 
	$counts->{nondup} + $counts->{duplicate};
exit; # bam files should automatically be closed





######################## Subroutines ###################################################

sub process_black_list {
	if ($black_list and -e $black_list) {
		eval {require Bio::ToolBox::Data};
		my $i = 0;
		eval {require Set::IntervalTree; $i = 1;};
		unless ($i) {
			warn " PROBLEM! Please install Set::IntervalTree to use black lists\n";
			undef $black_list;
			return;
		}
		my $Data = Bio::ToolBox::Data->new(file => $black_list) or 
			die "unable to read black list file '$black_list'\n";
		my %black_list_hash;
		$Data->iterate( sub {
			my $row = shift;
			$black_list_hash{ $row->seq_id } ||= [];
			push @{ $black_list_hash{ $row->seq_id } }, [$row->start - 1, $row->end + 0];
		} );
		printf " Loaded %d black list regions\n", ($Data->last_row);
		return \%black_list_hash;
	}
	return;
}



sub count_alignments {
	print " Counting $items";
	if ($parallel and $cpu > 1) {
		print " in $cpu threads....\n";
	}
	else {
		print "....\n";
	}
	
	# count the alignments in single or multi-thread
	my ($totalCount, $depth2count);
	if ($parallel and $cpu > 1) {
		($totalCount, $depth2count) = count_alignments_multithread();
	} else {
		($totalCount, $depth2count) = count_alignments_singlethread();
	}
	
	# report duplication statistics
		# remember that depth2count is simply depth and number of positions at that depth
		# not the number of alignments
		# to get total alignments sum( map {$_ * $depth2count{$_}} keys %depth2count )
	my @depths = sort {$a <=> $b} keys %$depth2count;
	my $nondupCount = sum(values %$depth2count);# assumes one unique at every single position
	my $dupCount = $totalCount - $nondupCount; # essentially count of positions with 2+ alignments
	my $dupRate = $dupCount / $totalCount;
	my $maxObserved = max(keys %$depth2count);
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
					$sumCount += ($depth2count->{$d} * $d);
				}
				else {
					# throw away the excess
					$sumCount += ($depth2count->{$d} * $depth);
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
			foreach my $d (keys %$depth2count) {
				$availableDups -= (($d - $max) * $depth2count->{$d}) 
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


sub count_alignments_singlethread {
	my $totalCount = 0;
	my %depth2count; # for generating a histogram of duplicate numbers
				   # key=depth, value=number of bases with this depth
	
	# count each chromosome one at a time
	for my $tid (@tid_list) {
	
		# prepare callback data structure
		my $data = {
			position    => -1, # a non-coordinate
			totalCount  => 0,  # total count of reads
			black_list  => undef,
			depth2count => {}, # depth to count hash
			reads       => [], # temp buffer for collecting reads
		};
	
		# process black lists for this chromosome
		# since we're using external interval tree module that is not fork-safe, must 
		# recreate interval tree each time
		my $seq_id = $sam->target_name($tid);
		if ($black_list_hash and exists $black_list_hash->{$seq_id}) {
			my $tree = Set::IntervalTree->new;
			foreach (@{ $black_list_hash->{$seq_id} }) {
				# don't need to insert any particular value, just want the interval
				$tree->insert(1, $_->[0], $_->[1]);
			}
			$data->{black_list} = $tree;
		}
		
		# walk through the reads on the chromosome
		my $seq_length = $sam->target_len($tid);
		$callback = $paired ? \&count_pe_callback : \&count_se_callback;
		low_level_bam_fetch($sam, $tid, 0, $seq_length, $callback, $data);
	
		# check to make sure we don't leave something behind
		if (defined $data->{reads}->[0]) {
			$paired ? count_up_pe_alignments($data) : count_up_se_alignments($data);
		}
		
		# add chromosome counts to global values
		$totalCount += $data->{totalCount};
		foreach my $d (keys %{$data->{depth2count}}) {
			$depth2count{$d} += $data->{depth2count}{$d};
		}
	}
	return ($totalCount, \%depth2count);
}


sub count_alignments_multithread {
	my $totalCount = 0;
	my %depth2count; # for generating a histogram of duplicate numbers
				   # key=depth, value=number of bases with this depth
	
	# set up manager
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_finish( sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
		# add child counts to global values
		$totalCount += $data->{totalCount};
		foreach my $d (keys %{$data->{depth2count}}) {
			$depth2count{$d} += $data->{depth2count}{$d};
		}
	});
	
	# run chromosomes in parallel
	for my $tid (@tid_list) {
		my $seq_length = $sam->target_len($tid);
		my $seq_id = $sam->target_name($tid);
		$pm->start and next;
		
		### in child
		$sam->clone;
		
		# prepare callback data structure
		my $data = {
			position    => -1, # a non-coordinate
			totalCount  => 0,  # total count of reads
			black_list  => undef,
			depth2count => {}, # depth to count hash
			reads       => [], # temp buffer for collecting reads
		};
	
		# process black lists for this chromosome
		# since we're using external interval tree module that is not fork-safe, must 
		# recreate interval tree each time
		if ($black_list_hash and exists $black_list_hash->{$seq_id}) {
			my $tree = Set::IntervalTree->new or die "unable to make interval tree!";
			foreach (@{ $black_list_hash->{$seq_id} }) {
				# don't need to insert any particular value, just want the interval
				$tree->insert(1, $_->[0], $_->[1]);
			}
			$data->{black_list} = $tree;
		}
		
		# walk through the reads on the chromosome
		$callback = $paired ? \&count_pe_callback : \&count_se_callback;
		low_level_bam_fetch($sam, $tid, 0, $seq_length, $callback, $data);
		
		# check to make sure we don't leave something behind
		if (defined $data->{reads}->[0]) {
			$paired ? count_up_pe_alignments($data) : count_up_se_alignments($data);
		}
		
		# finish and return to parent
		delete $data->{reads};
		delete $data->{position};
		undef $data->{black_list}; # why is this necessary?
		delete $data->{black_list};
		$pm->finish(0, $data); 
	}
	$pm->wait_all_children;
	return ($totalCount, \%depth2count);
}



### single-end alignment callback for reading
sub count_se_callback {
	my ($a, $data) = @_;
	
	# filter black listed regions
	if (defined $data->{black_list}) {
		my $results = $data->{black_list}->fetch($a->pos, $a->calend);
		return if @$results;
	}
	
	$data->{totalCount}++;
	
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
	
	# filter reads
	return unless $a->proper_pair; # consider only proper pair alignments
	return unless $a->tid == $a->mtid;
	return if $a->reversed; # count only forward alignments
	return unless $a->mreversed;
	return if $a->isize < 0; # wierd RF pair, a F alignment should only have + isize
	if (defined $data->{black_list}) {
		# filter black listed regions
		my $results = $data->{black_list}->fetch($a->pos, $a->calend);
		return if @$results;
	}
	
	$data->{totalCount}++;
	
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
		$data->{depth2count}{scalar @{ $fends{$pos} }} += 1;
	}
	foreach my $pos (keys %rends) {
		$data->{depth2count}{scalar @{ $rends{$pos} }} += 1;
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
		$data->{depth2count}{scalar @{$sizes{$s}}} += 1;
	}
}



### Deduplicate the bam file
sub deduplicate {
	
	print " Removing duplicates and writing new bam file";
	if ($parallel and $cpu > 1) {
		print " in $cpu threads....\n";
	}
	else {
		print "....\n";
	}
	
	# set callbacks and subroutines
	$callback = $paired ? \&write_pe_callback : \&write_se_callback;
	if ($random and defined $chance) {
		$write_out_alignments = $paired ? \&write_out_random_pe_alignments : 
			\&write_out_random_se_alignments;
	}
	else {
		$write_out_alignments = $paired ? \&write_out_max_pe_alignments : 
			\&write_out_max_se_alignments;
	}
	
	# deduplicate in single or multi-thread
	if ($parallel and $cpu > 1) {
		return deduplicate_multithread();
	} else {
		return deduplicate_singlethread();
	}
}


sub deduplicate_singlethread {
	
	# set up counters
	my %counts = (
		total      => 0,
		nondup     => 0,
		duplicate  => 0,
		toss       => 0,
		unpaired   => 0,
	);
	
	# open bam new bam file
	my $outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
		die "unable to open output bam file $outfile! $!";
		# using an unexported subroutine imported as necessary depending on bam availability
# 	my $outheader = $header->new; # make a new header object
# 	$outheader->text($htext); # add updated header text
# 	$outbam->header_write($outheader);
	$outbam->header_write($header);
	
	# set the random seed
	if (defined $seed) {
		srand($seed);
	}
	
	# walk through bam file again
	for my $tid (@tid_list) {
	
		# prepare callback data structure
		my $data = {
			total      => 0,
			nondup     => 0,
			duplicate  => 0,
			toss       => 0,
			position   => -1,
			outbam     => $outbam,
			black_list => undef,
			reads      => [],
			keepers    => {}, # hash of names of rev pe reads to keep
			dupkeepers => {}, # hash of names of rev pe dup reads to keep
		};
	
		# process black lists for this chromosome
		# since we're using external interval tree module that is not fork-safe, must 
		# recreate interval tree each time
		my $seq_id = $sam->target_name($tid);
		if ($black_list_hash and exists $black_list_hash->{$seq_id}) {
			my $tree = Set::IntervalTree->new;
			foreach (@{ $black_list_hash->{$seq_id} }) {
				# don't need to insert any particular value, just want the interval
				$tree->insert(1, $_->[0], $_->[1]);
			}
			$data->{black_list} = $tree;
		}
		
		# walk through the reads on the chromosome
		my $seq_length = $sam->target_len($tid);
		low_level_bam_fetch($sam, $tid, 0, $seq_length, $callback, $data);
	
		# check to make sure we don't leave something behind
		# possible with single-end, should be unlikely with paired-end
		if (defined $data->{reads}->[0]) {
			&$write_out_alignments($data);
		}
		if ($paired and scalar keys %{$data->{keepers}}) {
			$counts{unpaired} += scalar(keys %{$data->{keepers}});
		}
		
		# add up the counts
		printf " wrote %d alignments\n", $data->{nondup} + $data->{duplicate};
		foreach my $k (qw(total nondup duplicate toss)) {
			$counts{$k} += $data->{$k};
		}
	}
	
	return (\%counts, $outbam);
}


sub deduplicate_multithread {
	
	# set up counters
	my %counts = (
		total      => 0,
		nondup     => 0,
		duplicate  => 0,
		toss       => 0,
		unpaired   => 0,
	);
	
	# set up manager
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_finish( sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
		# add child counts to global values
		foreach my $k (qw(total nondup duplicate toss unpaired)) {
			$counts{$k} += $data->{$k};
		}
	});
	
	# generate child temporary bam files
	my $tempfile = $outfile;
	$tempfile =~ s/\.bam$//i;
	my %targetfiles = map {$_ => $tempfile . ".temp.$_.bam"} @tid_list;
	
	# walk through bam file again
	for my $tid (@tid_list) {
		$pm->start and next;
		
		### in child
		$sam->clone;
		
		# set the random seed
		if (defined $seed) {
			srand($seed);
		}
		
		# prepare a temporary bam file
		my $tf = $targetfiles{$tid};
		my $tempbam = Bio::ToolBox::db_helper::write_new_bam_file($tf) or 
			die "unable to open output bam file $tf! $!";
		$tempbam->header_write($header);
		
		# prepare callback data structure
		my $data = {
			total      => 0,
			nondup     => 0,
			duplicate  => 0,
			toss       => 0,
			unpaired   => 0,
			position   => -1,
			outbam     => $tempbam,
			black_list => undef,
			reads      => [],
			keepers    => {}, # hash of names of rev pe reads to keep
			dupkeepers => {}, # hash of names of rev pe dup reads to keep
		};
	
		# process black lists for this chromosome
		# since we're using external interval tree module that is not fork-safe, must 
		# recreate interval tree each time
		my $seq_id = $sam->target_name($tid);
		if ($black_list_hash and exists $black_list_hash->{$seq_id}) {
			my $tree = Set::IntervalTree->new;
			foreach (@{ $black_list_hash->{$seq_id} }) {
				# don't need to insert any particular value, just want the interval
				$tree->insert(1, $_->[0], $_->[1]);
			}
			$data->{black_list} = $tree;
		}
		
		# walk through the reads on the chromosome
		my $seq_length = $sam->target_len($tid);
		low_level_bam_fetch($sam, $tid, 0, $seq_length, $callback, $data);
	
		# check to make sure we don't leave something behind
		# possible with single-end, should be unlikely with paired-end
		if (defined $data->{reads}->[0]) {
			&$write_out_alignments($data);
		}
		if ($paired and scalar keys %{$data->{keepers}}) {
			$data->{unpaired} += scalar(keys %{$data->{keepers}});
		}
		
		# finish and return to parent
		delete $data->{outbam};
		delete $data->{reads};
		delete $data->{keepers};
		delete $data->{dupkeepers};
		delete $data->{position};
		undef $data->{black_list}; # why is this necessary?
		delete $data->{black_list};
		undef $tempbam;
		$pm->finish(0, $data); 
	}
	$pm->wait_all_children;
	
	
	# attempt to use external samtools to merge the bam files
	# this should be faster than going through Perl and bam adapters
	my $sam_app = which('samtools');
	if ($sam_app and not $no_sam) {
		my $command = sprintf "%s cat -o %s ", $sam_app, $outfile;
		$command .= join(' ', map { $targetfiles{$_} } @tid_list);
		print " executing $sam_app cat to merge children...\n";
		if (system($command)) {
			die "something went wrong with command '$command'!\n";
		}
		unlink values(%targetfiles);
		return (\%counts, $outbam);
	}
	
	
	# open final bam new bam file
	my $outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
		die "unable to open output bam file $outfile! $!";
		# using an unexported subroutine imported as necessary depending on bam availability
# 	my $href = ref($header);
# 	my $outheader = $href->new; # make a new header object
# 	$outheader->text($htext); # add updated header text
# 	$outbam->header_write($outheader);
# 	$header->text($htext); # add updated header text
	$outbam->header_write($header);
	
	# now remerge all the child files and write to main bam file
	if ($BAM_ADAPTER eq 'sam') {
		# open the low level bam file
		# avoid the usual high level bam file opening because that will force an 
		# automatic index, which we don't want or need
		# We will read the alignments as is. It should already be sorted and in order.
		# Must read the header first!
		for my $tid (@tid_list) {
			my $tf = $targetfiles{$tid};
			my $inbam = Bio::DB::Bam->open($tf) or die "unable to open $tf!\n";
			my $h = $inbam->header;
			while (my $a = $inbam->read1) {
				&$write_alignment($outbam, $a);
			}
			undef $inbam;
			unlink $tf;
		} 
	}	
	elsif ($BAM_ADAPTER eq 'hts') {
		for my $tid (@tid_list) {
			my $tf = $targetfiles{$tid};
			my $inbam = Bio::DB::HTSfile->open($tf) or die "unable to open $tf!\n";
			my $h = $inbam->header_read;
			while (my $a = $inbam->read1($h)) {
				&$write_alignment($outbam, $a);
			}
			undef $inbam;
			unlink $tf;
		}
	}
	
	# finished
	return (\%counts, $outbam);
}

### single-end alignment callback for writing
sub write_se_callback {
	my ($a, $data) = @_;
	
	# filter black listed regions
	if (defined $data->{black_list}) {
		my $results = $data->{black_list}->fetch($a->pos, $a->calend);
		return if @$results;
	}
	
	$data->{total}++;
	
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
	if (defined $data->{black_list}) {
		# filter black listed regions
		my $results = $data->{black_list}->fetch($a->pos, $a->calend);
		return if @$results;
	}
	
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
	$data->{total}++ unless $a->reversed; # only count forward alignments
	
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
		$data->{nondup}++;
		
		# always write one alignment
		my $a1 = shift @{ $fends{$pos} };
		&$write_alignment($data->{outbam}, $a1);
		
		# write remaining alignments up to max
		for (my $i = 1; $i < $max; $i++) {
			my $a = shift @{ $fends{$pos} };
			last unless $a;
			&$write_alignment($data->{outbam}, $a);
			$data->{duplicate}++;
		}
		$data->{toss} += scalar @{ $fends{$pos} };
		if ($mark and scalar @{ $fends{$pos} }) {
			while (my $a = shift @{ $fends{$pos} }) {
				mark_alignment($a);
				&$write_alignment($data->{outbam}, $a);
			}
		}
	}
	
	# write reverse reads
	foreach my $pos (keys %rends) {
		$data->{nondup}++;
		
		# always write one alignment
		my $a1 = shift @{ $rends{$pos} };
		&$write_alignment($data->{outbam}, $a1);
		
		# write remaining alignments up to max
		for (my $i = 1; $i < $max; $i++) {
			my $a = shift @{ $rends{$pos} };
			last unless $a;
			&$write_alignment($data->{outbam}, $a);
			$data->{duplicate}++;
		}
		$data->{toss} += scalar @{ $rends{$pos} };
		if ($mark and scalar @{ $rends{$pos} }) {
			while (my $a = shift @{ $rends{$pos} }) {
				mark_alignment($a);
				&$write_alignment($data->{outbam}, $a);
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
		$data->{nondup}++;
		
		# always write one alignment
		my $a1 = shift @{ $f_sizes{$s} };
		&$write_alignment($data->{outbam}, $a1);
		$data->{keepers}{$a1->qname} = 1; # remember name for reverse read
		
		# write remaining up to max
		for (my $i = 1; $i < $max; $i++) {
			my $a = shift @{ $f_sizes{$s} };
			last unless $a;
			&$write_alignment($data->{outbam}, $a);
			$data->{keepers}{$a->qname} = 1; # remember name for reverse read
			$data->{duplicate}++;
		}
		$data->{toss} += scalar @{ $f_sizes{$s} };
		if ($mark and scalar @{ $f_sizes{$s} }) {
			while (my $a = shift @{ $f_sizes{$s} }) {
				mark_alignment($a);
				&$write_alignment($data->{outbam}, $a);
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
				&$write_alignment($data->{outbam}, $a);
				delete $data->{keepers}{$name};
			}
			elsif (exists $data->{dupkeepers}{$name}) {
				# we have processed the forward alignment as a duplicate keeper
				# immediately write the reverse alignment
				&$write_alignment($data->{outbam}, $a);
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
			$data->{nondup}++;
			&$write_alignment($data->{outbam}, $fends{$pos}->[0]);
		}
		else {
			random_se_alignment_writer($data, $fends{$pos});
		}
	}
	
	# write reverse reads
	foreach my $pos (keys %rends) {
		if (scalar @{ $rends{$pos} } == 1) {
			# there's only one alignment, automatic keep
			$data->{nondup}++;
			&$write_alignment($data->{outbam}, $rends{$pos}->[0]);
		}
		else {
			random_se_alignment_writer($data, $rends{$pos});
		}
	}

}


sub random_se_alignment_writer {
	my ($data, $alignments) = @_;
	
	# always write the first one
	$data->{nondup}++;
	&$write_alignment($data->{outbam}, shift @{$alignments});
	
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
	$data->{duplicate} += scalar(@keep);
	foreach my $a (@keep) {
		&$write_alignment($data->{outbam}, $a);
	}
	
	# record losers
	$data->{toss} += scalar(@toss);
	if ($mark) {
		foreach my $a (@toss) {
			mark_alignment($a);
			&$write_alignment($data->{outbam}, $a);
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
			$data->{nondup}++;
			my $a = $f_sizes{$s}->[0];
			&$write_alignment($data->{outbam}, $a);
			$data->{keepers}{$a->qname} = 1; # remember name for reverse read
		}
		else {
			# always write the first one
			$data->{nondup}++;
			my $a = shift @{ $f_sizes{$s} };
			&$write_alignment($data->{outbam}, $a);
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
			$data->{duplicate} += scalar(@keep);
			foreach my $a (@keep) {
				&$write_alignment($data->{outbam}, $a);
				$data->{keepers}{$a->qname} = 1; # remember name to write reverse read
			}
	
			# record losers
			$data->{toss} += scalar(@toss);
			if ($mark) {
				foreach my $a (@toss) {
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
					$data->{dupkeepers}{$a->qname} = 1; # remember name 
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
				&$write_alignment($data->{outbam}, $a);
				delete $data->{keepers}{$name};
			}
			elsif (exists $data->{dupkeepers}{$name}) {
				# we have processed the forward alignment as a duplicate keeper
				# immediately write the reverse alignment
				&$write_alignment($data->{outbam}, $a);
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

sub mark_alignment {
	# mark alignments as a duplicate
	my $f = $_[0]->flag;
	unless ($f & 0x400) {
		$f += 0x400;
		$_[0]->flag($f);
	}
}



