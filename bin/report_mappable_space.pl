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
use List::Util qw(sum);
use Set::IntSpan::Fast;
use Parallel::ForkManager;
use Bio::ToolBox::db_helper 1.50 qw(
	open_db_connection
	low_level_bam_fetch
	$BAM_ADAPTER
);

our $VERSION = 1.2;

my @bamfiles;
my $min_mapq = 10;
my $cpu      = 4;
my $chr_exclude;

unless (@ARGV) {
	print <<END;

A script to report the mappable space of a genome based on empirical mapping 
derived from one or more bam files. For a given ChIPSeq experiment, provide all 
available Bam files.

Two numbers are reported: all mappable space reported by all alignments, 
regardless of mapping quality and number of hits to the genome, and unique 
mappability, as determined by a minimum mapping quality score (default $min_mapq).
Results are written to standard out.

USAGE:  report_mappable_space.pl *.bam
        report_mappable_space.pl --chrskip '^chrRandom.+|chrM' *.bam

VERSION: $VERSION
       
OPTIONS:
  -i --in <file>         An input bam file, must be sorted and indexed
                           Repeat for each input file
  -q --qual <int>        Minimum mapping quality to be considered unique ($min_mapq)
  --chrskip <regex>      Provide a regular expression for skipping unwanted chromosomes
  -c --cpu <int>         Specify the number of threads to use ($cpu)

END
	exit;
}

GetOptions(
	'i|in=s'        => \@bamfiles,       # the input bam file path
	'q|qual|mapq=i' => \$min_mapq,       # minimum mapping quality
	'chrskip=s'     => \$chr_exclude,    # skip chromosomes
	'c|cpu=i'       => \$cpu,            # number of cpu cores to use
	'bam=s'         => \$BAM_ADAPTER,    # specifically set the bam adapter, advanced!
) or die "unrecognized parameters\n";

if ( scalar(@ARGV) ) {
	push @bamfiles, @ARGV;
}
unless (@bamfiles) {
	die " At least one bam file must be provided!\n";
}
foreach my $f (@bamfiles) {
	unless ( -e $f and -r _ ) {
		die " Input Bam file '$f' does not exist or can't be read!\n";
	}
	unless ( $f =~ /\.bam$/ ) {
		die " File '$f' not a Bam file! Only Bam files are supported\n";
	}
}

# global values
my $unique_space = 0;
my $all_space    = 0;
my %chr_list;

### Generate chromosome list
# use first bam file as an example
my $example = $bamfiles[0];
{
	my $sam = open_db_connection( $example, 1 )
		or die "unable to read bam file '$example'!\n";
	for my $tid ( 0 .. $sam->n_targets - 1 ) {
		my $chr = $sam->target_name($tid);
		my $len = $sam->target_len($tid);
		if ( $chr_exclude and $chr =~ /$chr_exclude/xi ) {
			next;
		}
		$chr_list{$chr} = $len;
	}
	undef $sam;
}

### Initialize fork manager
my $pm = Parallel::ForkManager->new($cpu)
	or die "unable to initialize ForkManager object!\n";
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) = @_;
		$unique_space += $data->{unique};
		$all_space    += $data->{all};
	}
);

### Calculate mappability space
foreach my $chr ( keys %chr_list ) {
	$pm->start and next;

	# initialize interval counters
	my $unique_set = Set::IntSpan::Fast->new()
		or die "unable to initialize Integer Set!\n";
	my $all_set = Set::IntSpan::Fast->new();
	my $data    = {
		unique => $unique_set,
		all    => $all_set
	};

	# loop through each bam file
	my $length = $chr_list{$chr};
	foreach my $file (@bamfiles) {

		# make sure we do not use a cached database object!!!
		my $sam = open_db_connection( $file, 1 )
			or die "unable to open bam file '$file'!\n";

		# check chromosome and length
		my ( $tid, undef, undef ) = $sam->header->parse_region($chr);
		my $seq_length = $sam->target_len($tid);
		if ( $length != $seq_length ) {
			warn
" ! Sequence length discrepancy for $chr in $file! Compare $length with $seq_length\n";
		}

		# re-initialize coordinate checker
		$data->{start} = 0;
		$data->{end}   = 0;

		# scan the chromosome
		low_level_bam_fetch( $sam, $tid, 0, $seq_length, \&callback, $data );

		# finished with this file
		undef $sam;
	}

	# calculate space for unique
	my $unique_sum = 0;
	my $iter       = $unique_set->iterate_runs();
	while ( my ( $from, $to ) = $iter->() ) {
		$unique_sum += ( $to - $from + 1 );
	}

	# calculate space for all
	my $all_sum = 0;
	$iter = $all_set->iterate_runs();
	while ( my ( $from, $to ) = $iter->() ) {
		$all_sum += ( $to - $from + 1 );
	}

	# finished with this chromosome
	$data->{unique} = $unique_sum;
	$data->{all}    = $all_sum;
	$pm->finish( 0, $data );
}
$pm->wait_all_children;

### Finish
my $genome = sum( values %chr_list );
printf
"\n Total Genome: %.3f Mb\n All mappable space: %.3f Mb (%d%%)\n Unique mappable space %.3f Mb (%d%%)\n",
	( $genome / 1000000 ), ( $all_space / 1000000 ), ( ( $all_space / $genome ) * 100 ),
	( $unique_space / 1000000 ), ( ( $unique_space / $genome ) * 100 );

### Alignment callback
sub callback {
	my ( $a, $data ) = @_;

	# coordinates
	my $s = $a->pos + 1;
	my $e = $a->calend;

	# discard duplicates on the assumption that they don't add anything
	return if ( $s == $data->{start} and $e == $data->{end} );
	$data->{start} = $s;
	$data->{end}   = $e;

	# check if this range has been covered
	# presumption that test is easier and faster than adding new range
	# particularly for subsequent files and the set has a full chromosome already
	# there's a minuscule chance it might not be in unique but is in all
	return if $data->{all}->contains_all_range( $s, $e );

	# record
	$data->{all}->add_range( $s, $e );
	if ( $a->qual >= $min_mapq ) {
		$data->{unique}->add_range( $s, $e );
	}
}

