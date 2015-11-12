#!/usr/bin/env perl

use strict;
use Bio::ToolBox::db_helper qw(get_chromosome_list);

unless (@ARGV) {
	print <<USAGE;
A script to print out the lengths of chromosomes present in a database.
A database can be Bio::DB::SeqFeature::Store database, Bam file (.bam),
bigWig (.bw) file, bigBed (.bb) file, fasta (.fa or .fasta) file, or 
a directory of individual fasta files. 

Usage: $0 <database>

Chromosome names and sizes in bp are printed to standard out.

USAGE
	exit;
}


# print chromosome name and lengths
my $database = shift @ARGV;
my @chromosomes = get_chromosome_list($database);
foreach (@chromosomes) {
	print "$_->[0]\t$_->[1]\n";
}
