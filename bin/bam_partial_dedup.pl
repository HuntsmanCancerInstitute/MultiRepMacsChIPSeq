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
use File::Which;
use IO::File;
use List::Util qw(sum min max);
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

my $VERSION = 5.2;

my $DOC = <<END;

A script to remove excessive duplicate alignments down to an acceptable 
fraction. This is in contrast to traditional duplicate removers that 
remove all duplicates, retaining only one alignment per position. 

Duplicate reads may be artificial (derived from PCR duplication due to 
very low input) or natural (multiple DNA fragments enriched at few narrow 
foci, such as a ChIP peaks). Removing all duplicates can significantly 
reduce high-enrichment peaks, but not removing duplicates can lead to 
false positives. An optimal balance is therefore desirable. This is 
most important when comparing between replicates or samples.

This script can randomly subsample or remove duplicate reads to reach a 
target duplication rate. This results in a more uniform duplicate reduction 
across the genome, which can be more typical of true PCR duplication. 
Set the target duplication fraction rate with the --frac option below. 

This script can also simply remove excessive duplicate reads at positions 
that exceed a specified target threshold. This can be set either alone or
in combination with the random subsample. Using alone is generally not 
recommmended, as it reduces signal at extreme peaks without addressing 
low level duplication elsewhere across the genome.

For ChIPSeq applications, check the duplication level of the input. For 
mammalian genomes, typically 5-20% duplication is observed in sonicated 
input. For very strong enrichment of certain targets, it's not unusual 
to see higher duplication rates in ChIP samples than in Input samples. 
Generally, set the target fraction of all samples to the lowest observed 
duplication rate.

Single-end aligment duplicates are checked for start position, strand, and 
calculated alignment end position to check for duplicates. Because of this, 
the numbers may be slightly different than calculated by traditional duplicate 
removers.

Paired-end alignments are treated as fragments. Only properly paired 
alignments are considered; singletons are skipped and dropped. Fragments 
are checked for start position and fragment length (paired insertion size) 
for duplicates. Random subsampling should not result in broken pairs.

Optical duplicates, arising from neighboring clusters on a sequencing flow 
cell with identical sequence, may now be checked. When random subsampling 
duplicates, optical duplicates should critically be ignored. This is highly 
recommended for patterned flow cells from Illumina NovaSeq or NextSeq. Set a 
distance of 100 pixels for unpatterned (Illumina HiSeq) or at least 2500 for  
patterned (NovaSeq). By default, optical duplicate alignments are not written 
to output. To ONLY filter for optical duplicates, set --max to a very high number.
Note that tile-edge duplicates are not counted as such.

Bam files representing multiple sequence lanes that are identified properly by 
unique read groups are de-deduplicated separately with regard to optical 
duplicates, but not coordinate duplicates. (New as of version 3)

Existing alignment duplicate marks (bit flag 0x400) are ignored. 

Since repetitive and high copy genomic regions are a big source of duplicate 
alignments, these regions can and should be entirely skipped by providing a 
file with recognizable coordinates. Any alignments overlapping these intervals 
are skipped in both counting and writing. 

USAGE:  bam_partial_dedup.pl --in input.bam
        bam_partial_dedup.pl --frac 0.xx --in input.bam --out output.bam

VERSION: $VERSION
       
OPTIONS:
  --in <file>         The input bam file, should be sorted and indexed
  --out <file>        The output bam file containing unique and retained 
                        duplicates; optional if you're just checking the 
                        duplication rate.
  --pe                Bam files contain paired-end alignments and only 
                        properly paired duplicate fragments will be checked for 
                        duplication. Singletons are silently dropped.
  --qual <int>        Skip alignments below indicated mapping quality (0)
  --mark              Write non-optical duplicate alignments to output marked 
                        with flag bit 0x400.
  --frac <float>      Decimal fraction representing the target duplication 
                        rate in the final file. 
  --max <int>         Integer representing the maximum number of alignments 
                        at each position. Set to 1 to remove all duplicates.
  --optical           Enable optical duplicate checking
  --distance <int>    Set optical duplicate distance threshold.
                        Use 100 for unpatterned flowcell (HiSeq) or 
                        2500 for patterned flowcell (NovaSeq). Default 100.
                        Setting this value automatically sets --optical.
  --report            Write duplicate distance report files only, no de-duplication
  --keepoptical       Keep optical duplicates in output as marked 
                        duplicates with flag bit 0x400. Optical duplicates 
                        are not differentiated from non-optical duplicates.
  --coord <string>    Provide the tile:X:Y integer 1-base positions in the 
                        read name for optical checking. For Illumina CASAVA 1.8 
                        7-element names, this is 5:6:7 (default)
  --blacklist <file>  Provide a bed/gff/text file of repeat regions to skip
  --chrskip <regex>   Provide a regex for skipping certain chromosomes
  --seed <int>        Provide an integer to set the random seed generator to 
                        make the subsampling consistent (non-random).
  --cpu <int>         Specify the number of threads to use (4) 
  --verbose           Print more information
  --help              Print full documentation

END

unless (@ARGV) {
	print <<END;

 A script for checking and completely or partially removing duplicate 
 alignments from a Bam file.
 
 VERSION: $VERSION
 
 SIMPLE USAGE:
 Check single-end duplication rate:
   bam_partial_dedup.pl --in input.bam
 
 Check paired-end duplication rate:
   bam_partial_dedup.pl --pe --in input.bam
 
 Random subsample single-end duplicates to 10%:
   bam_partial_dedup.pl --frac 0.1 --in input.bam --out output.bam

 Use -h or --help for full documentation and options.
 
END
	exit;
}

### Get Options
my ($fraction, $max, $infile, $outfile, $do_optical, $optical_thresh, $keep_optical, 
	$mark, $random, $paired, $min_mapq, $chr_exclude, $black_list, $seed, $cpu, $tilepos, 
	$xpos, $ypos, $name_coordinates, $report_distance, $no_sam, $verbose, $help);
my @program_options = @ARGV;
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'o|out=s'    => \$outfile, # name of output file 
	'optical!'   => \$do_optical, # check for optical duplicates
	'distance=i' => \$optical_thresh, # optical threshold distance
	'report!'    => \$report_distance, # write duplicate distance report
	'keepoptical!' => \$keep_optical, # include optical duplicates in output
	'mark!'      => \$mark, # mark the duplicates
	'random!'    => \$random, # flag to random downsample duplicates - obsolete
	'frac=f'     => \$fraction, # target fraction of duplicates
	'max=i'      => \$max, # the maximum number of alignments per position
	'pe!'        => \$paired, # treat as paired-end alignments
	'qual|mapq=i' => \$min_mapq, # minimum mapping quality
	'chrskip=s'  => \$chr_exclude, # skip chromosomes
	'blacklist=s' => \$black_list, # file of high copy repeat regions to avoid
	'seed=i'     => \$seed, # seed for non-random random subsampling
	'cpu=i'      => \$cpu, # number of cpu cores to use
	'coord=s'    => \$name_coordinates, # tile:X:Y name positions
	'verbose!'   => \$verbose, # tell me more
	'bam=s'      => \$BAM_ADAPTER, # specifically set the bam adapter, advanced!
	'nosam!'     => \$no_sam, # avoid using external sam adapter, advanced!
	'help!'      => \$help, # print documentation
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

my $max_optical_rate = 0.01;

### Check options
if ($help) {
	print $DOC;
	exit 0;
}
unless ($infile) {
	die " No input file provided!\n";
}
if (defined $max) {
	if ($max < 1) {
		die " Maximum number of alignments must be at least 1!\n";
	}
	elsif ($max == 1) {
		# no need to calculate fractions if user set max number to 1
		undef $fraction;
	}
}
if ($fraction and $fraction > 1) {
	die " Unrecognized fraction '$fraction'. Should be decimal (0.x)\n";
}
if ($random) {
	print " Specifying option --random is no longer necessary and is default\n";
}
if ($report_distance) {
	# optical must be on
	print " Enabling optical duplicate checking when reporting\n" if not $do_optical;
	$do_optical = 1;
}
if ($optical_thresh) {
	$do_optical = 1;
}
if ($do_optical) {
	# default for HiSeq non-patterned flow cells, perhaps we should use 2500 for NovaSeq
	unless ($optical_thresh) {
		print " Setting optical distance default to 100 pixels\n";
		$optical_thresh = 100;
	}
}
if ($name_coordinates) {
	if ($name_coordinates =~ /(\d):(\d):(\d)/) {
		# coordinates must be converted to 0-based indices
		$tilepos = $1 - 1;
		$xpos = $2 - 1;
		$ypos = $3 - 1;
	}
	else {
		die " Name coordinate must be integers as tile:X:Y, such as 5:6:7\n";
	}
}
else {
	# these defaults are for Illumina CASAVA 1.8+ Fastq data in 0-based coordinates
	$tilepos = 4;
	$xpos = 5;
	$ypos = 6;
}
if ($parallel and not defined $cpu) {
	$cpu = 4;
}
if ($min_mapq) {
	unless ($min_mapq > 0 and $min_mapq < 256) {
		die " Invalid minimum mapping quality! Must be an integer between 0-255!\n";
	}
}


### Open bam files
# input bam file
my $sam = open_db_connection($infile) or # automatically takes care of indexing
	die "unable to read bam $infile $!";
print " using '$BAM_ADAPTER' Perl Bam adapter\n" if $verbose;

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
my $chance = 0; # probability for subsampling
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
exit 0 unless $outfile;
exit 0 if $report_distance;
exit 0 unless (defined $max or $chance > 0);
$outfile .= '.bam' unless $outfile =~ /\.bam$/i;
# update header
my $htext = $header->text;
$htext .= sprintf("\@PG\tID:bam_partial_dedup\tVN:%s\tCL:%s", $VERSION, $0);
$htext .= " --pe" if $paired;
$htext .= " --max=$max" if ($max);
$htext .= " --frac=$fraction" if ($fraction);
$htext .= " --seed=$seed" if $seed;
$htext .= " --mark" if $mark;
$htext .= " --optical --distance=$optical_thresh " if ($do_optical);
$htext .= " --keepoptical" if $keep_optical;
$htext .= " --chrskip='$chr_exclude'" if $chr_exclude;
$htext .= " --blacklist=$black_list" if $black_list;
$htext .= " --qual=$min_mapq" if $min_mapq;
$htext .= " --coord $name_coordinates" if $name_coordinates;
$htext .= " --bam $BAM_ADAPTER --in $infile --out $outfile\n";

my ($counts, $outbam) = deduplicate();



### Print results and finish
printf "  Total mapped: %34d 
  Optical duplicate count: %23d
  Non-duplicate count: %27d
  Retained non-optical duplicate count: %10d
  Removed non-optical duplicate count: %11d
  Current duplication rate: %22.4f\n",
	$counts->{total}, $counts->{optical}, $counts->{nondup}, $counts->{duplicate}, 
	$counts->{toss}, $counts->{duplicate}/($counts->{duplicate} + $counts->{nondup});
if ($paired and $counts->{unpaired}) {
	printf "  Unpaired count: %32d\n", $counts->{unpaired};
}
printf " Wrote %d alignments to $outfile\n", $paired ? 
	($counts->{nondup} + $counts->{duplicate}) * 2 : 
	$counts->{nondup} + $counts->{duplicate};
exit 0; # bam files should automatically be closed





######################## Subroutines ###################################################

sub process_black_list {
	return unless $black_list;
	if (-e $black_list) {
		eval {require Bio::ToolBox::Data}; # this should be available
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
	else {
		print " PROBLEM: Specified black list file '$black_list' not found, skipping\n";
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
	my ($totalCount, $opticalCount, $coordErrorCount, $depth2count, $dupmatrix, $distance2count);
	if ($parallel and $cpu > 1) {
		($totalCount, $opticalCount, $coordErrorCount, $depth2count, $dupmatrix, 
			$distance2count) = count_alignments_multithread();
	} else {
		($totalCount, $opticalCount, $coordErrorCount, $depth2count, $dupmatrix,
			 $distance2count) = count_alignments_singlethread();
	}
	
	# sanity check
	if ($totalCount == 0) {
		print "\n Something is wrong! Zero $items were counted!\n Check your bam file to make sure mapped alignments are present.\n";
		if ($paired) {
			print " Also check that the bam file truly has paired-end alignments,\n or run in single-end mode (without --pe option)\n\n";
		}
		exit;
	}
	
	# report duplication statistics
		# remember that depth2count is simply depth and number of positions at that depth
		# not the number of alignments
		# to get total alignments sum( map {$_ * $depth2count{$_}} keys %depth2count )
	my @depths = sort {$a <=> $b} keys %$depth2count;
	my $opticalRate = $opticalCount / $totalCount;
	my $workingCount = $totalCount - $opticalCount; # non-optical dup count
	my $nondupCount = sum(values %$depth2count);# assumes one unique at every single position
	my $dupCount = $workingCount - $nondupCount; # essentially count of positions with 2+ alignments
	my $dupRate = $dupCount / $workingCount;
	my $maxObserved = max(keys %$depth2count);
	# right justify the numbers in the printf to make it look pretty
	printf "  Total mapped: %24d
  Optical duplicate count: %13d
  Optical duplicate rate: %14.4f
  Non-optical working count: %11d
  Non-optical duplicate count: %9d
  Non-optical duplication rate: %8.4f
  Non-duplicate count: %17d
  Maximum position depth: %14d
  Mean position depth: %17.4f\n", 
		$totalCount, $opticalCount, $opticalRate, $workingCount, 
		$dupCount, $dupRate, $nondupCount, $maxObserved, $workingCount / $nondupCount;
	
	# warning about optical duplicate errors
	if ($coordErrorCount) {
		printf "  WARNING: %d Alignments did not have read name pixel coordinates for optical detection!\n",
			$coordErrorCount;
		if ( $coordErrorCount > int($totalCount/10) ) {
			# this count is greater than 10%, optical de-duplication not reliable
			print "   Error count is > 10%, turning off optical deduplication\n";
			$do_optical = 0;
		}
	}
	
	# write distance reports
	if ($report_distance) {
		# first get the limits of the matrix
		my @xvalues;
		my @yvalues;
		while (my ($x, $yhash) = each %$dupmatrix) {
			push @xvalues, $x;
			push @yvalues, max(keys %$yhash);
		}
		my $xlimit = min(max(@xvalues), 50000);
		my $ylimit = min(max(@yvalues), 50000);
		
		# write matrix to file
		my $matrixfile = $outfile || $infile;
		$matrixfile =~ s/\.bam$//;
		$matrixfile .= '.dup-distance.matrix.txt';
		my $fh = IO::File->new($matrixfile, '>') or 
			die "unable to write to file $matrixfile! $!";
		$fh->print("# duplicate delta pixel distances for $infile\n");
		$fh->print("# number is duplicate count observed at X and Y delta distance from previous\n");
		$fh->printf("Y\t%s\n", join("\t", map {$_ * 50} (0..$xlimit)));
		foreach my $i (0..$ylimit) {
			$fh->printf("%d\t%s\n", $i * 50, 
				join("\t", map { $dupmatrix->{$_}{$i} || 0 } (0..$xlimit) ) 
			);
		}
		$fh->close;
		print "  wrote duplicate delta distance matrix to file $matrixfile\n";
		
		# write max distance hash to file
		my $histogramfile = $outfile || $infile;
		$histogramfile =~ s/\.bam$//;
		$histogramfile .= '.dup-threshold.txt';
		my $fh = IO::File->new($histogramfile, '>') or 
			die "unable to write to file $histogramfile! $!";
		$fh->print("# Fraction of observed duplicates at the given pixel distance threshold\n");
		my $dup_sum = $opticalCount + $dupCount;
		$fh->print("# Total observed duplicates for $infile: $dup_sum\n");
		my $current_dup_sum = $dup_sum;
		$fh->print("Distance\t$infile\n");
		for my $i (0..400) {
			$fh->printf("%d\t%0.04f\n", $i * 50, ($current_dup_sum/$dup_sum));
			$current_dup_sum -= $distance2count->{$i} || 0;
		}
		$fh->close;
		print "  wrote table of duplicate fraction at each pixel distance to file $histogramfile\n";
	}
	
	# check if we need to continue
	if ($fraction and $dupRate <= $fraction) {
		# the duplication rate is acceptable, no de-duplication necessary
		print " Actual duplication rate is less than indicated target fraction ($fraction).\n";
		
		# check to see whether max or optical duplication is still necessary
		if (defined $max and $max >= 1 and $do_optical) {
			# maximum depth is set plus checking for opticals
			
			if ($maxObserved <= $max) {
				print " Maximum observed depth is less than maximum allowed ($max).\n";
				if ($opticalRate > $max_optical_rate) {
					print " Optical duplicate rate is higher than allowed ($max_optical_rate).\n";
					# proceed with the given max value, de-duplication will proceed
				}
				else {
					# both optical rate and max duplicate is OK
					# reset the max allowed value, de-duplication will not proceed
					print " No de-duplication necessary.\n";
					undef $max;
				}
			}
			else {
				# we will continue with max cutoff
				print " Maximum observed depth is above maximum allowed ($max).\n";
			}
		}
		elsif (defined $max and $max >= 1 and not $do_optical) {
			if ($maxObserved <= $max) {
				print " Maximum observed depth is less than maximum allowed ($max).\n No de-duplication necessary.\n";
				# both optical rate and max duplicate is OK
				# reset the max allowed value, de-duplication will not proceed
				undef $max;
			}
			else {
				# we will continue with max cutoff
				print " Maximum observed depth is above maximum allowed ($max).\n";
			}
		
		}
		elsif (not defined $max and $do_optical) {
			# no maximum was set, duplication is acceptable
			# check optical rate
			if ($opticalRate > $max_optical_rate) {
				print " Optical duplicate rate is higher than allowed ($max_optical_rate).\n";
				# de-duplication will proceed
				# artificially set the max value to force the de-duplication routine
				# without actually removing non-optical duplicates
				$max = $maxObserved + 1;
			}
			else {
				print " Optical duplication rate is ok.\n No de-duplication necessary.\n";
				# optical rate is OK
				# de-duplication will not proceed
			}
		}
	}
	
	# check maximum
	elsif (defined $max and $maxObserved <= $max and not $fraction) {
		print " Maximum observed alignment depth is less than maximum allowed ($max).\n";
		# check for optical rate
		if ($do_optical and $opticalRate > $max_optical_rate) {
			print " Optical duplicate rate is higher than allowed ($max_optical_rate).\n";
			# de-duplication will proceed
			# artificially set the max value to force the de-duplication routine
			# without actually removing non-optical duplicates
			$max = $maxObserved + 1;
		}
		else {
			# no de-duplication is necessary
			print " No de-duplication necessary\n";
			undef $max;
		}
	}
	
	# calculate fraction
	elsif ($fraction) {
		print " Calculating expectations for random subsampling\n";
		# this is the number of duplicates available to be subsampled
		my $availableDups = $dupCount;
		if (defined $max) {
			# adjust for the number of duplicates to be tossed for exceeding max cutoff
			while (my ($d, $c) = each %$depth2count) {
				if ($max > 1 and $d > $max) {
					$availableDups -= (($d - $max) * $c) 
				}
			}
		}
		
		# this is what the new final total should be
		my $newTotal = int($nondupCount / (1 - $fraction) );
		printf "  Expected new total: %18d\n", $newTotal;
		
		# this is how many can be duplicates
		my $allowedDups = int($newTotal * $fraction);
		printf "  Expected duplicates kept: %12d\n", $allowedDups;
		
		# this is the chance for being kept
		$chance = sprintf("%.8f", $allowedDups / $availableDups); 
		printf "  Probability of keeping: %14s\n", $chance;
	}
	
	# otherwise we simply return
}


sub count_alignments_singlethread {
	my $totalCount = 0;
	my $opticalCount = 0;
	my $coordErrorCount = 0;
	my %depth2count; # for generating a histogram of duplicate numbers
				   # key=depth, value=number of bases with this depth
	my %dupmatrix; # xy matrix for duplicate reporting
	my %distance2count; # single axis distance histogram reporting
	
	# count each chromosome one at a time
	for my $tid (@tid_list) {
	
		# prepare callback data structure
		my $data = {
			position       => -1, # a non-coordinate
			totalCount     => 0,  # total count of reads
			optCount       => 0, # optical count
			coordError     => 0, # pixel coordinate errors
			black_list     => undef,
			depth2count    => {}, # depth to count hash
			dupmatrix      => {}, # duplicate x,y coordinates for reporting
			distance2count => {}, # max distance to count hash
			reads          => [], # temp buffer for collecting reads
		};
	
		# process black lists for this chromosome
		# since we're using external interval tree module that is not fork-safe, must 
		# recreate interval tree each time
		my $seq_id = $sam->target_name($tid);
		print "  counting $seq_id $items...\n" if $verbose;
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
		$opticalCount += $data->{optCount};
		$coordErrorCount += $data->{coordError};
		printf("   %d $items counted on $seq_id\n", $data->{totalCount}) if $verbose;
		while (my ($d, $c) = each %{ $data->{depth2count} } ) {
			$depth2count{$d} += $c;
		}
		
		# duplicate distances
		if ($report_distance) {
			# duplicate x,y matrix
			foreach my $x (keys %{ $data->{dupmatrix} }) {
				$dupmatrix{$x} ||= {};
				foreach my $y ( keys %{ $data->{dupmatrix}{$x} } ) {
					$dupmatrix{$x}{$y} += $data->{dupmatrix}{$x}{$y};
				}
			}
			# distance to key hash
			while (my ($d, $c) = each %{ $data->{distance2count} }) {
				$distance2count{$d} += $c;
			}
		}
		
		# pixel coordinate errors
		$coordErrorCount += $data->{coordError};
	}
	return ($totalCount, $opticalCount, $coordErrorCount, \%depth2count, 
		\%dupmatrix, \%distance2count);
}


sub count_alignments_multithread {
	my $totalCount = 0;
	my $opticalCount = 0;
	my $coordErrorCount = 0;
	my %depth2count; # for generating a histogram of duplicate numbers
				   # key=depth, value=number of bases with this depth
	my %dupmatrix; # xy matrix for duplicate reporting
	my %distance2count; # single axis distance histogram reporting
	
	# set up manager
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_finish( sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
		# add child counts to global values
		print "   child process $pid exited with code $exit_code\n" if $verbose;
		$totalCount += $data->{totalCount};
		$opticalCount += $data->{optCount};
		$coordErrorCount += $data->{coordError};
		while (my ($d, $c) = each %{ $data->{depth2count} } ) {
			$depth2count{$d} += $c;
		}
		if ($report_distance) {
			# duplicate x,y matrix
			foreach my $x (keys %{ $data->{dupmatrix} }) {
				$dupmatrix{$x} ||= {};
				foreach my $y ( keys %{ $data->{dupmatrix}{$x} } ) {
					$dupmatrix{$x}{$y} += $data->{dupmatrix}{$x}{$y};
				}
			}
			# distance to key hash
			while (my ($d, $c) = each %{ $data->{distance2count} }) {
				$distance2count{$d} += $c;
			}
		}
	});
	
	# run chromosomes in parallel
	for my $tid (@tid_list) {
		my $seq_length = $sam->target_len($tid);
		my $seq_id = $sam->target_name($tid);
		$pm->start and next;
		
		### in child
		$sam->clone;
		print "  counting $seq_id $items...\n" if $verbose;
		
		# prepare callback data structure
		my $data = {
			position       => -1, # a non-coordinate
			totalCount     => 0,  # total count of reads
			optCount       => 0, # optical count
			coordError     => 0, # pixel coordinate errors
			black_list     => undef,
			depth2count    => {}, # depth to count hash
			dupmatrix      => {}, # duplicate x,y coordinates for reporting
			distance2count => {}, # max distance to count hash
			reads          => [], # temp buffer for collecting reads
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
		printf("   %d $items counted on $seq_id\n", $data->{totalCount}) if $verbose;
		delete $data->{reads};
		delete $data->{position};
		undef $data->{black_list}; # why is this necessary?
		delete $data->{black_list};
		$pm->finish(0, $data); 
	}
	$pm->wait_all_children;
	
	return ($totalCount, $opticalCount, $coordErrorCount, \%depth2count, 
		\%dupmatrix, \%distance2count);
}



### single-end alignment callback for reading
sub count_se_callback {
	my ($a, $data) = @_;
	
	# mapping quality
	return if ($min_mapq and $a->qual < $min_mapq); # mapping quality
	
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
	return if ($min_mapq and $a->qual < $min_mapq); # mapping quality
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
		# count up based on optical test
		if ($do_optical) {
			my ($nondup, $dup, $coord_error) = 
				identify_optical_duplicates($fends{$pos}, $data);
			$data->{depth2count}{scalar @$nondup} += 1;
			$data->{optCount} += scalar @$dup;
			$data->{coordError} += $coord_error;
		}
		else {
			# blindly take everything
			$data->{depth2count}{ scalar @{$fends{$pos}} } += 1;
		}
	}
	foreach my $pos (keys %rends) {
		# count up based on optical test
		if ($do_optical) {
			my ($nondup, $dup, $coord_error) = 
				identify_optical_duplicates($rends{$pos}, $data);
			$data->{depth2count}{scalar @$nondup} += 1;
			$data->{optCount} += scalar @$dup;
			$data->{coordError} += $coord_error;
		}
		else {
			# blindly take everything
			$data->{depth2count}{ scalar @{$rends{$pos}} } += 1;
		}
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
		# count up based on optical test
		if ($do_optical) {
			my ($nondup, $dup, $coord_error) = identify_optical_duplicates($sizes{$s}, $data);
			$data->{depth2count}{scalar @$nondup} += 1;
			$data->{optCount} += scalar @$dup;
			$data->{coordError} += $coord_error;
		}
		else {
			# blindly take everything
			$data->{depth2count}{ scalar @{$sizes{$s}} } += 1;
		}
	}
}



### Identify optical duplicates
sub identify_optical_duplicates {
	my $alignments = shift;
	my $data = shift || undef; # this is only for reporting distances
	
	# check for only one alignment and quickly return
	if (scalar @$alignments == 1) {
		return ($alignments, []);
	}
	
	# extract tile and coordinates from each alignment given and sort 
	my %tile2alignment;
	my @dups;
	my @nondups;
	my $coord_error = 0;
	foreach my $a (@$alignments) {
		my @bits = split(':', $a->qname);
		if (
			$#bits < $ypos or (
				not defined $bits[$tilepos] and 
				not defined $bits[$xpos] and 
				not defined $bits[$ypos]
			)
		) {
			# not able to identify coordinates in the alignment name
			# record error count and place in nondups array
			$coord_error++;
			push @nondups, $a;
		}
		else {
			my $rg = $a->aux_get("RG") || '1';
			my $rgtile = join('.', $rg, $bits[$tilepos]);
			$tile2alignment{$rgtile} ||= [];
			# each item in the array is another array of [x, y, $a]
			push @{$tile2alignment{$rgtile}}, [$bits[$xpos], $bits[$ypos], $a];
		}
	}
	
	# check whether we have actionable alignments
	if (scalar keys %tile2alignment == 0) {
		# couldn't find any coordinates!
		return (\@nondups, \@dups, $coord_error);
	}
	
	
	
	# walk through each read-group tile
	while (my ($rgtile, $rgtile_alnts) = each %tile2alignment) {
		# check the number of alignments
		if (scalar @$rgtile_alnts == 1) {
			# we only have one, so it's a nondup
			push @nondups, $rgtile_alnts->[0][2];
		}
		else {
			# we have more than one so must sort through
			
			# sort by increasing X coordinate
			my @spots = sort {$a->[0] <=> $b->[0]} @$rgtile_alnts;
			# collect the deltas from one spot to the next on the X axis
			# start with second sorted spot and subtract the one before it
			my @xdiffs = map { $spots[$_]->[0] - $spots[$_-1]->[0] } (1..$#spots);
			
			if ($report_distance) {
				# we need to generate a report, so take the extra time to calculate the
				# absolute y deltas using the same sorted order of spots
				my @ydiffs = map { abs($spots[$_]->[1] - $spots[$_-1]->[1]) } (1..$#spots);
				for my $i (0..$#xdiffs) {
					# note that we're recording based on xdiffs array, not spots array
					# therefore we're not recording the first spot in the spots array
					# i.e. one non-duplicate will be left out
					my $x = int($xdiffs[$i] / 50);
					my $y = int($ydiffs[$i] / 50);
					$data->{dupmatrix}{$x} ||= {};
					$data->{dupmatrix}{$x}{$y} += 1;
					
					# record at what pixel cutoff this would be considered an optical duplicate
					my $m = max($x, $y);
					# an optical duplicate would be within this distance in both axes
					# looking at up to 20,000 pixels is reasonable range to check
					if ($m <= 400) {
						$data->{distance2count}{$m} += 1;
					}
					
				}
			}
			
			
			## Check 
			# working on the assumption here that by definition an optical cluster is 
			# going to be below the threshold on both axes, and not just one axis
			# therefore we only need to check the x axis first
			if (min(@xdiffs) > $optical_thresh) {
				# everything is ok
				foreach (@spots) {
					push @nondups, $_->[2];
				}
			}
			else {
				# two or more spots are too close on X axis, must check out
				my $first = shift @spots;
				while (@spots) {
					if ($spots[0]->[0] - $first->[0] < $optical_thresh) {
						# second is too close to first
						# keep taking spots until they're no longer close
						my @closest;
						push @closest, $first, shift @spots;
						# continue comparing the next one with the previous one
						while (
							scalar @spots and 
							$spots[0]->[0] - $closest[-1]->[0] < $optical_thresh
						) {
							push @closest, shift @spots;
						}
						
						## Now must process this X cluster
						# check Y coordinates by sorting and calculating deltas
						@closest = sort {$a->[1] <=> $b->[1]} @closest;
						my @ydiffs = map { $closest[$_]->[1] - $closest[$_-1]->[1] } (1..$#closest);
						
						# check the Y deltas in this X cluster
						if (min(@ydiffs) > $optical_thresh) {
							# they're all good
							foreach (@closest) {
								push @nondups, $_->[2];
							}
						}
						else {
							# we definitely have some XY clusters within threshold
							my $xyfirst = shift @closest;
							while (@closest) {
								if ($closest[0]->[1] - $xyfirst->[1] < $optical_thresh) {
									# these two are close
									my @clustered;
									push @clustered, $xyfirst, shift @closest;
									# continue compareing the next one with the previous one
									while (
										scalar @closest and 
										$closest[0]->[1] - $clustered[-1]->[1] < $optical_thresh
									) {
										push @clustered, shift @closest;
									}
									# done, no more
									
									# take the first one as pseudo random chosen one
									push @nondups, $clustered[0]->[2];
									# remainder are optical duplicates
									for my $i (1..$#clustered) {
										push @dups, $clustered[$i]->[2];
									}
									
									## prepare for next round
									$xyfirst = shift @closest || undef;
								}
								else {
									# this one is ok
									push @nondups, $xyfirst->[2];
									$xyfirst = shift @closest;
								}
							}
							
							# check for last remaining alignment
							push @nondups, $xyfirst->[2] if defined $xyfirst;
						}
						
						## Prepare for next round
						$first = shift(@spots) || undef;
					}
					else {
						# first is ok
						push @nondups, $first->[2];
						$first = shift @spots;
					}
					# continue
				}
				
				# check for last remaining alignment
				push @nondups, $first->[2] if defined $first;
			}
		}
	}
	
	# finished
	return (\@nondups, \@dups, $coord_error);
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
	if (defined $chance) {
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
		optical    => 0,
		nondup     => 0,
		duplicate  => 0,
		toss       => 0,
		unpaired   => 0,
	);
	
	# open bam new bam file
	my $outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
		die "unable to open output bam file $outfile! $!";
		# using an unexported subroutine imported as necessary depending on bam availability
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
			optical    => 0,
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
		print "  deduplicating $seq_id $items...\n" if $verbose;
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
		printf("   wrote %d $items for $seq_id\n", $data->{nondup} + $data->{duplicate}) 
			if $verbose;
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
		optical    => 0,
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
		print "   child process $pid exited with code $exit_code\n" if $verbose;
		foreach my $k (qw(total optical nondup duplicate toss unpaired)) {
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
		my $seq_id = $sam->target_name($tid);
		print "  deduplicating $seq_id $items...\n" if $verbose;
		
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
			optical    => 0,
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
		printf("   wrote %d $items for $seq_id\n", $data->{nondup} + $data->{duplicate}) 
			if $verbose;
		$pm->finish(0, $data); 
	}
	$pm->wait_all_children;
	
	
	# attempt to use external samtools to merge the bam files
	# this should be faster than going through Perl and bam adapters
	my $sam_app = which('samtools');
	if ($sam_app and not $no_sam) {
		
		# write temporary new header file
		my $samfile = $outfile;
		$samfile =~ s/\.bam$/temp.sam/;
		my $fh = IO::File->new($samfile, '>') or 
			die "unable to write temporary sam file\n";
		$fh->print($htext);
		$fh->close;
		
		# check samtools version
		my $sam_check = qx($sam_app --version);
		my $sam_version;
		if ($sam_check =~ /samtools 1\.(\d+)/) {
			$sam_version = $1;
		}
		
		# run external samtools concatenate
		my $command = sprintf("%s cat -h %s -o %s ", $sam_app, $samfile, $outfile);
		if ($sam_version >= 10) {
			# enable multi-theading and no program note for 1.10 or greater
			$command .= sprintf("--no-PG --threads %s ", $cpu);
		}
		$command .= join(' ', map { $targetfiles{$_} } @tid_list);
		print " executing $sam_app cat to merge children...\n";
		if (system($command)) {
			die "something went wrong with command '$command'!\n";
		}
		unlink $samfile, values(%targetfiles);
		return (\%counts, $outbam);
	}
	
	
	# open final bam new bam file
	print " writing merged Bam file with adapter\n";
	my $outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
		die "unable to open output bam file $outfile! $!";
		# using an unexported subroutine imported as necessary depending on bam availability
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
	
	# mapping quality
	return if ($min_mapq and $a->qual < $min_mapq); # mapping quality
	
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
	return if ($min_mapq and $a->qual < $min_mapq); # mapping quality
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
	while (my ($pos, $ends) = each %fends) {
		
		# collect the non-optical duplicate reads
		my $reads;
		if ($do_optical) {
			my ($nondup, $dup, $coord_error) = identify_optical_duplicates($ends);
			$reads = $nondup;
			$data->{optical} += scalar @$dup;
			if ($keep_optical and scalar(@$dup)) {
				# go ahead and write out optical duplicates now
				foreach my $a (@$dup) {
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
				}
			}
		}
		else {
			# blindly take everything
			$reads = $ends;
		}
		
		# always write one alignment
		$data->{nondup}++;
		my $a1 = shift @$reads;
		&$write_alignment($data->{outbam}, $a1);
		
		# write remaining alignments up to max
		for (my $i = 1; $i < $max; $i++) {
			my $a = shift @$reads;
			last unless $a;
			&$write_alignment($data->{outbam}, $a);
			$data->{duplicate}++;
		}
		$data->{toss} += scalar @$reads;
		if ($mark and scalar @$reads) {
			while (my $a = shift @$reads) {
				mark_alignment($a);
				&$write_alignment($data->{outbam}, $a);
			}
		}
	}
	
	# write reverse reads
	while (my ($pos, $ends) = each %rends) {
		
		# collect the non-optical duplicate reads
		my $reads;
		if ($do_optical) {
			my ($nondup, $dup, $coord_error) = identify_optical_duplicates($ends);
			$reads = $nondup;
			$data->{optical} += scalar @$dup;
			if ($keep_optical and scalar(@$dup)) {
				# go ahead and write out optical duplicates now
				foreach my $a (@$dup) {
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
				}
			}
		}
		else {
			# blindly take everything
			$reads = $ends;
		}
		
		# always write one alignment
		$data->{nondup}++;
		my $a1 = shift @$reads;
		&$write_alignment($data->{outbam}, $a1);
		
		# write remaining alignments up to max
		for (my $i = 1; $i < $max; $i++) {
			my $a = shift @$reads;
			last unless $a;
			&$write_alignment($data->{outbam}, $a);
			$data->{duplicate}++;
		}
		$data->{toss} += scalar @$reads;
		if ($mark and scalar @$reads) {
			while (my $a = shift @$reads) {
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
		
		# collect the non-optical duplicate reads
		my @reads;
		if ($do_optical) {
			my ($nondup, $dup, $coord_error) = identify_optical_duplicates($f_sizes{$s});
			@reads = @$nondup;
			$data->{optical} += scalar @$dup;
			if ($keep_optical and scalar(@$dup)) {
				# go ahead and write out optical duplicates now
				foreach my $a (@$dup) {
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
					$data->{dupkeepers}{$a->qname} = 1; # remember name for rev
				}
			}
		}
		else {
			# blindly take everything
			@reads = @{$f_sizes{$s}};
		}
		
		# always write one alignment
		$data->{nondup}++;
		my $a1 = shift @reads;
		&$write_alignment($data->{outbam}, $a1);
		$data->{keepers}{$a1->qname} = 1; # remember name for reverse read
		
		# write remaining up to max
		for (my $i = 1; $i < $max; $i++) {
			my $a = shift @reads;
			last unless $a;
			&$write_alignment($data->{outbam}, $a);
			$data->{keepers}{$a->qname} = 1; # remember name for reverse read
			$data->{duplicate}++;
		}
		$data->{toss} += scalar @reads;
		if ($mark and scalar @reads) {
			while (my $a = shift @reads) {
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
				mark_alignment($a);
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
		
		# collect the non-optical duplicate reads
		my @reads;
		if ($do_optical) {
			my ($nondup, $dup, $coord_error) = identify_optical_duplicates($fends{$pos});
			@reads = @$nondup;
			$data->{optical} += scalar @$dup;
			if ($keep_optical and scalar(@$dup)) {
				# go ahead and write out optical duplicates now
				foreach my $a (@$dup) {
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
				}
			}
		}
		else {
			# blindly take everything
			@reads = @{$fends{$pos}};
		}
		
		
		if (scalar @reads == 1) {
			# there's only one alignment, automatic keep
			$data->{nondup}++;
			&$write_alignment($data->{outbam}, $reads[0]);
		}
		else {
			random_se_alignment_writer($data, \@reads);
		}
	}
	
	# write reverse reads
	foreach my $pos (keys %rends) {
		
		# collect the non-optical duplicate reads
		my @reads;
		if ($do_optical) {
			my ($nondup, $dup, $coord_error) = identify_optical_duplicates($rends{$pos});
			@reads = @$nondup;
			$data->{optical} += scalar @$dup;
			if ($keep_optical and scalar(@$dup)) {
				# go ahead and write out optical duplicates now
				foreach my $a (@$dup) {
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
				}
			}
		}
		else {
			# blindly take everything
			@reads = @{$rends{$pos}};
		}
		
		# write out reads
		if (scalar @reads == 1) {
			# there's only one alignment, automatic keep
			$data->{nondup}++;
			&$write_alignment($data->{outbam}, $reads[0]);
		}
		else {
			random_se_alignment_writer($data, \@reads);
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
	if (defined $max and scalar(@keep) >= $max) {
		push @toss, splice(@keep, $max - 1); # we've already written the first one
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
		
		# collect the non-optical duplicate reads
		my @reads;
		if ($do_optical) {
			my ($nondup, $dup, $coord_error) = identify_optical_duplicates($f_sizes{$s});
			@reads = @$nondup;
			$data->{optical} += scalar @$dup;
			if ($keep_optical and scalar(@$dup)) {
				# go ahead and write out optical duplicates now
				foreach my $a (@$dup) {
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
					$data->{dupkeepers}{$a->qname} = 1; # remember name for rev
				}
			}
		}
		else {
			# blindly take everything
			@reads = @{$f_sizes{$s}};
		}
		
		# write out reads
		if (scalar @reads == 1) {
			# there's only one fragment, automatic keep
			$data->{nondup}++;
			&$write_alignment($data->{outbam}, $reads[0]);
			$data->{keepers}{$reads[0]->qname} = 1; # remember name for reverse read
		}
		else {
			# always write the first one
			$data->{nondup}++;
			my $a = shift @reads;
			&$write_alignment($data->{outbam}, $a);
			$data->{keepers}{$a->qname} = 1; # remember name for reverse read
			
			# duplicates get to roll the dice
			my @keep;
			my @toss;
			while (my $a = shift @reads) {
				if (rand(1) < $chance) {
					push @keep, $a;
				}
				else {
					push @toss, $a;
				}
			}
	
			# cut off above maximum
			if (defined $max and scalar(@keep) >= $max) {
				push @toss, splice(@keep, $max - 1); # we've already written the first one
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
				mark_alignment($a);
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



