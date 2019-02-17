#!/usr/bin/env perl

use strict;
use File::Copy;
use Getopt::Long;
use Bio::ToolBox::big_helper qw(generate_chromosome_file);

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
  -o --out "file"              Optional file name

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

$db ||= shift @ARGV;
unless ($db) {
	die "must specify a database file!\n";
}


# generate the chromosome file
# this is always a randomly named temp file
my $chromo_file = generate_chromosome_file($db, $chrskip) or 
	die "unable to generate chromosome file!\n";

# if requested, rename to specified output file
if ($out) {
	move($chromo_file, $out);
	print "wrote chromosome file '$out'\n";
}
else {
	print "wrote chromosome file '$chromo_file'\n";
}

