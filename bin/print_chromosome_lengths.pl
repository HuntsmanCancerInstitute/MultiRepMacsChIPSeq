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
use English qw(-no_match_vars);
use File::Spec;
use IO::File;
use Getopt::Long;
use Bio::ToolBox::db_helper 1.69 qw(get_chromosome_list);

our $VERSION = 1.1;

unless (@ARGV) {
	print <<USAGE;
A script to print out the lengths of chromosomes present in a database.
A database can be Bio::DB::SeqFeature::Store database, Bam file (.bam),
bigWig (.bw) file, bigBed (.bb) file, multi-fasta (.fa or .fasta) file, or 
a directory of individual fasta files. 

If more than one source file is given, then all are checked for consistency 
of chromosome names and order. An output file is only written if source 
files have the same chromosome list.

Version: $VERSION

Usage: print_chromosome_lengths.pl <database1> ...
       print_chromosome_lengths.pl -K 'chrM|alt|contig|un' -o chrom.sizes *.bam

Options:
  -d --db "file"               Indexed database file, may repeat
  -K --chrskip "text"          Chromosome skip regex
  -o --out "file"              Output file name, default STDOUT

USAGE
	exit;
}

### Options
my @databases;
my $chrskip;
my $out;
GetOptions(
	'd|db=s'      => \@databases,
	'K|chrskip=s' => \$chrskip,
	'o|out=s'     => \$out
) or die "unrecognized option!\n";

# database file
if ( scalar(@ARGV) ) {
	push @databases, @ARGV;
}
unless (@databases) {
	die "must specify one or more database files!\n";
}

### Collect Chromosomes
my %seqstring;
my @chrlist;

foreach my $db (@databases) {

	# get the chromosomes for this database
	my @chromosomes = get_chromosome_list( $db, $chrskip );
	unless (@chromosomes) {
		warn " no chromosome sequences in database file '$db'!\n";
		next;
	}

	# check chromosome string
	my $string = join( "\t", map { $_->[0] } @chromosomes );
	$seqstring{$string} ||= [];
	push @{ $seqstring{$string} }, $db;

	# final chromosome list
	unless (@chrlist) {
		@chrlist = @chromosomes;
	}
}

### Validate chromosomes
if ( scalar( keys %seqstring ) > 1 ) {
	printf
"\n WARNING!!! There are %d different chromosome name and/or order lists\n present in the source files!!! This may/will present problems with\n downstream applications! This may be correctable without re-aligning.\n",
		scalar( keys %seqstring );
	foreach my $k ( keys %seqstring ) {
		printf "\n=== Files:\n  %s\n=== Chromosomes:\n  %s\n",
			join( "\n  ", @{ $seqstring{$k} } ),
			join( "\n  ", split( /\t/, $k ) );
	}
	exit 0;
}

# write the chromosome list
if ($out) {
	my $fh = IO::File->new( $out, '>' )
		or die "unable to write file! $OS_ERROR\n";
	foreach my $chr (@chrlist) {
		$fh->printf( "%s\t%d\n", $chr->[0], $chr->[1] );
	}
	$fh->close;
	print " Successfully wrote chromosome file '$out'\n";
}
else {
	foreach my $chr (@chrlist) {
		printf( "%s\t%d\n", $chr->[0], $chr->[1] );
	}
}

