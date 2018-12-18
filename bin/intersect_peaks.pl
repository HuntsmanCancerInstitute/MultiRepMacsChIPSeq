#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Which;
use Bio::ToolBox::Data;

# variables
my $tool = which('bedtools');
my $outfile;
my @files;
my $help;


### Documentation
my $docs = <<DOC;

A script to intersect two or more peak files. This is a wrapper around the bedtools 
program.

It will first merge all of the peak files into a single representative bed file.

It will then run the bedtools Jaccard statistic pairwise across all of the peak 
files and write out a merged table of the results. The Jaccard statistic measures 
the amount of spatial overlap between two peak files (intersection/union) reported 
as a fraction between 0 and 1.

Three files will be written:
    basename.bed                    the merged peaks
    basename.jaccard.txt            the Jaccard results in a table
    basename.n_intersection.txt     the number of intersections in a table

USAGE: intersect_peaks.pl --out <basename> peak1.narrowPeak peak2.narrowPeak ....

OPTIONS:
    --out basename          Provide the output basename
    --bed path              Path to bedtools ($tool)
    --help                  Print documentation
DOC




### Options
GetOptions(
	'out=s'             => \$outfile,
	'bed=s'             => \$tool,
	'help!'             => \$help,
) or die "unrecognized option!\n";

# check options
if ($help) {
	print $docs;
	exit;
}
unless (@ARGV) {
	print $docs;
	exit;
}
die "must provide output base name!\n" unless $outfile;
die "bedtools not in your PATH!\n" unless $tool;
$outfile =~ s/\.(?:bed|narrowPeak)$//i; # strip any existing extension if provided

# inputs
my @files = @ARGV;
die "must provide 2 or more files!\n" unless scalar(@files) > 1;
my @names = @files;
foreach (@names) {s/\.narrowPeak//} 




### Merge the peaks
my $command = "cat ";
foreach my $f (@files) {
	$command .= sprintf("%s ", $f);
}
$command .= sprintf(" | %s sort -i - | %s merge -i - > %s",
	$tool, $tool, $outfile . '.bed' );
if (system($command)) {
	die "something went wrong! command:\n $command\n";
}



### Calculate Jaccard statistic
my $JaccardData = Bio::ToolBox::Data->new(
	columns => ['File', @names],
);
my $IntersectionData = Bio::ToolBox::Data->new(
	columns => ['File', @names],
);

for my $f1 (0..$#names) {
	my $j = $JaccardData->add_row();
	my $i = $IntersectionData->add_row();
	
	$JaccardData->value($j, 0, $names[$f1]);
	$IntersectionData->value($i, 0, $names[$f1]);
	
	for my $f2 (0..$#names) {
		my $result = `$tool jaccard -a $files[$f1] -b $files[$f2]`;
		my @lines = split("\n", $result);
		my (undef, undef, $jaccard, $number) = split(/\s+/, $lines[1]);
		$JaccardData->value($j, $f2 + 1, $jaccard);
		$IntersectionData->value($i, $f2 + 1, $number);
	}
}

$JaccardData->save("$outfile\.jaccard.txt");
$IntersectionData->save("$outfile\.n_intersection.txt");

