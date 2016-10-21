#!/usr/bin/perl

use strict;
use List::Util qw(sum max);
use Bio::ToolBox::db_helper::bam; 

unless (@ARGV) {
	print <<END;

A script to remove excessive duplicate alignments, leaving up to X number 
of alignments at any given position, where X can be assigned by the user. 
This is in contrast to traditional duplicate removers that retain only 
one alignment per position.

Currently only works with single-end data, for now. Paired-end alignments are 
treated like single end, likely breaking pairs. Alignments are not selected 
for retention; only the first X are retained.

The limit can be either an integer, where up to X alignments are tolerated,
or it can be a fraction (0.xx), whereupon the maximum number of alignments 
that achieve this target fraction is calculated.

NOTE: This de-duplication is slightly more aggressive than samtools rmdup, 
since it only examines the alignment start position and does not take into 
account variations in the end position due to clipping and read trimming.

Usage: $0 <limit> <input.bam> <output.bam>

END
	exit;
}

# max limit of reads at any given position
my ($target_fraction, $max);
my $limit = shift @ARGV;
if ($limit =~ /^\d+$/) {
	$max = $limit;
}
elsif ($limit =~ /^0\.\d+$/) {
	$target_fraction = $limit;
}
else {
	die "unrecognized limit '$limit'. Should be integer (X) or decimal (0.x)\n";
}



### Open bam files
# input bam file
my $infile = shift @ARGV;
my $sam = open_bam_db($infile) or # automatically takes care of indexing
	die "unable to read bam $infile $!";

# output bam file
my $outfile = shift @ARGV or die "no output file provided!\n";

# write header
my $header = $sam->bam->header;



# counters
my $totalCount = 0;
my %depth2count; # for generating a histogram of duplicate numbers
               # key=depth, value=number of bases with this depth


### count the reads and number at each position
if ($target_fraction) {
	# count alignments on each chromosome
	print " Counting alignments....\n";
	for my $tid (0 .. $sam->n_targets - 1) {
		my $seq_length = $sam->target_len($tid);
	
		# prepare callback data structure
		my $data = {
			position   => -1,
			reads      => [],
		};
	
		# walk through the reads on the chromosome
		$sam->bam_index->fetch($sam->bam, $tid, 0, $seq_length, \&count_callback, $data);
	
		# check to make sure we don't leave something behind
		if (defined $data->{reads}->[0]) {
			my $number = scalar @{ $data->{reads} };
			$depth2count{$number} += 1;
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
		if ($rate < $target_fraction) {
			$max = $depths[$i]; 
			last;
		}
	}
	printf "
 %12s total mapped reads
 %12s predicted duplicate reads
 %12s non-duplicate reads
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

### write out new bam file with specified number of targets
print " Removing duplicates and writing new bam file....\n";

# open bam new bam file
my $outbam = write_new_bam_file($outfile) or 
	die "unable to open output bam file $outfile! $!";
	# this uses low level Bio::DB::Bam object
$outbam->header_write($header);

$totalCount    = 0;
my $keepCount  = 0;
my $tossCount  = 0;
for my $tid (0 .. $sam->n_targets - 1) {
	my $seq_length = $sam->target_len($tid);
	
	# prepare callback data structure
	my $data = {
		position   => -1,
		reads      => [],
	};
	
	# walk through the reads on the chromosome
	$sam->bam_index->fetch($sam->bam, $tid, 0, $seq_length, \&write_callback, $data);
	
	# check to make sure we don't leave something behind
	if ($data->{reads}->[0]) {
		for (my $i = 0; $i < $limit; $i++) {
			my $a = shift @{ $data->{reads} };
			last unless $a;
			$outbam->write1($a);
		}
	}
}

# print results
printf "
 %12s total mapped reads
 %12s duplicate mapped reads discarded
 %12s mapped reads retained
 new duplicate rate %.4f
", $totalCount, $keepCount, $tossCount, $tossCount/$totalCount;

# finish up
undef $outbam;
check_bam_index($outfile);


### alignment callback for reading
sub count_callback {
	my ($a, $data) = @_;
	$totalCount++;
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		if ($data->{position} != -1) {
			my $number = scalar @{ $data->{reads} };
			$depth2count{$number} += 1;
		}
		$data->{position} = $pos;
		$data->{reads} = []; # clear out existing alignments
		push @{ $data->{reads} }, $a;
	}
}

### alignment callback for writing
sub write_callback {
	my ($a, $data) = @_;
	$totalCount++;
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		# write out the max number of reads
		for (my $i = 0; $i < $max; $i++) {
			my $a = shift @{ $data->{reads} };
			last unless $a;
			$outbam->write1($a);
			$keepCount++;
		}
		$tossCount += scalar @{ $data->{reads} };
		
		# reset
		$data->{reads} = [];
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}


