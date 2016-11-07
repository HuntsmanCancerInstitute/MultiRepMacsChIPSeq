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

The limit can be either an integer, where up to X alignments are tolerated,
or it can be a fraction (0.xx), whereupon the maximum number of alignments 
that achieve this target fraction is calculated.

Duplicates are checked for start position, strand, and calculated alignment 
end position to check for duplicates. Because of this, the numbers are slightly 
different than calculated by traditional duplicate removers.

Currently only works with single-end data, for now. Paired-end alignments are 
treated like single end, likely breaking pairs. Alignments are not selected 
for retention; only the first X are retained.


USAGE: bam_partial_dedup.pl <limit> <input.bam> <output.bam> 
       bam_partial_dedup.pl <0.xx> <input.bam> 
       bam_partial_dedup.pl <limit> <input.bam> <output.bam> <duplicate.bam>

If you simply want to determine the cutoff for a given duplication fraction, 
provide a fractional limit and omit the output file name.

If you want to retain the duplicates, provide an optional third bam file.

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
my $outfile = shift @ARGV || undef;

my $duplicatefile = shift @ARGV || undef;

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
			position   => -1, # a non-coordinate
			reads      => [],
		};
	
		# walk through the reads on the chromosome
		$sam->bam_index->fetch($sam->bam, $tid, 0, $seq_length, \&count_callback, $data);
	
		# check to make sure we don't leave something behind
		if (defined $data->{reads}->[0]) {
			count_up_alignments($data);
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
if ($duplicatefile) {
	my $dupbam = write_new_bam_file($duplicatefile) or 
		die "unable to open duplicate bam file $duplicatefile! $!";
		# this uses low level Bio::DB::Bam object
	$dupbam->header_write($header);
}

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
	if (defined $data->{reads}->[0]) {
		write_out_alignments($data);
	}
}

# print results
printf "
 %12s total mapped reads
 %12s duplicate mapped reads discarded
 %12s mapped reads retained
 new duplicate rate %.4f
", $totalCount, $tossCount, $keepCount, $tossCount/$totalCount;

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
			# a real position
			count_up_alignments($data);
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
		write_out_alignments($data);
		
		# reset
		$data->{reads} = [];
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}


sub count_up_alignments {
	my $data = shift;
	
	# split up based on end point and strand
	my %fends; # forward ends
	my %rends; # reverse ends
	foreach my $r (@{$data->{reads}}) {
		my $end = $r->calend; 
		if ($r->reversed) {
			$rends{$end} ||= [];
			push @{ $rends{$end} }, $r;
		}
		else {
			$fends{$end} ||= [];
			push @{ $fends{$end} }, $r;
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


sub write_out_alignments {
	my $data = shift;
	
	# split up based on end point and strand
	my %fends; # forward ends
	my %rends; # reverse ends
	foreach my $r (@{$data->{reads}}) {
		my $end = $r->calend; 
		if ($r->reversed) {
			$rends{$end} ||= [];
			push @{ $rends{$end} }, $r;
		}
		else {
			$fends{$end} ||= [];
			push @{ $fends{$end} }, $r;
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
		if ($dupbam and @{ $fends{$pos} }) {
			while (@{ $fends{$pos} }) {
				$dupbam->write1(shift @{ $fends{$pos} } );
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
		if ($dupbam and @{ $rends{$pos} }) {
			while (@{ $rends{$pos} }) {
				$dupbam->write1(shift @{ $rends{$pos} } );
			}
		}
	}
}


