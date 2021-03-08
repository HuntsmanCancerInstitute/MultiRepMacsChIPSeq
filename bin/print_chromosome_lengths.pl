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
use File::Spec;
use IO::File;
use Getopt::Long;
use Bio::ToolBox::db_helper qw(get_chromosome_list);

unless (@ARGV) {
	print <<USAGE;
A script to print out the lengths of chromosomes present in a database.
A database can be Bio::DB::SeqFeature::Store database, Bam file (.bam),
bigWig (.bw) file, bigBed (.bb) file, fasta (.fa or .fasta) file, or 
a directory of individual fasta files. 

Usage: $0 <database>

Options:
  -d --db "file"               Indexed database file
  -K --chrskip "text"          Chromosome skip regex
  -o --out "file"              Output file name, default db basename

USAGE
	exit;
}


### Options
my $db;
my $chrskip;
my $out;
GetOptions(
	'd|db=s'      => \$db,
	'K|chrskip=s' => \$chrskip,
	'o|out=s'     => \$out
) or die "unrecognized option!\n";

# database file
$db ||= shift @ARGV;
unless ($db) {
	die "must specify a database file!\n";
}

# output file
unless ($out) {
	(undef, undef, $out) = File::Spec->splitpath($db);
	$out =~ s/\.(?:bam|bw|bigwig|bb|bigbed|fasta|fa)$//i;
	$out .= '.chromosomes.txt';
}


# collect the chromosomes
my @chromosomes = get_chromosome_list($db, $chrskip);
unless (@chromosomes) {
	die " no chromosome sequences identified in database file $db!\n";
}


# write the chromosome list
my $fh = IO::File->new($out, '>') or 
	die "unable to write file! $!";
foreach my $chr (@chromosomes) {
	$fh->printf("%s\t%d\n", $chr->[0], $chr->[1]);
}
$fh->close;
print "wrote chromosome file '$out'\n";


