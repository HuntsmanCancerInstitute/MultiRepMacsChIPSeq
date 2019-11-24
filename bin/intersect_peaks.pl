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
use File::Spec;
use File::Which;
use Bio::ToolBox 1.65;

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

Finally, it will run bedtools multiinter tool to perform a multi-way intersection 
and calculate the intervals for each category of overlap. This is parsed into a 
summary file suitable for drawing a Venn diagram based on spatial overlap.

Five files will be written:
    basename.bed                    the merged peaks
    basename.jaccard.txt            the Jaccard results in a table
    basename.n_intersection.txt     the number of intersections in a table
    basename.multi.txt              data file from multi-intersection 
    basename.spatialVenn.txt        summary of spatial overlap for each category

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
my @names;
foreach my $f (@files) {
	# strip path and extension
	my (undef, $dir, $name) = File::Spec->splitpath($f);
	$name =~ s/\.(?:bed|narrowpeak|gappedpeak|broadpeak)(?:\.gz)?$//i;
	push @names, $name;
}




### Merge the peaks
print " Merging peak files....\n";
my $command = "cat ";
foreach my $f (@files) {
	$command .= sprintf("%s ", $f);
}
$command .= sprintf(" | %s sort -i - | %s merge -i -c 4,5 -o distinct,mean - > %s",
	$tool, $tool, $outfile . '.bed' );
if (system($command)) {
	die "something went wrong! command:\n $command\n";
}



### Calculate Jaccard statistic
print " Calculating Jaccard overlap....\n";
my $JaccardData = Bio::ToolBox->new_data('File', @names);
my $IntersectionData = Bio::ToolBox->new_data('File', @names);

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




### Multiple intersection
print " Calculating multi-way intersection....\n";
my $multi_file = $outfile . '.multi.txt';
$command = sprintf("%s multiinter -header -names %s -i %s > %s", 
	$tool, join(' ', @names), join(' ', @files), $multi_file);
if (system($command)) {
	die "something went wrong! command:\n $command\n";
}

# parse lengths
my $total = 0;
my %name2space;
my $MultiData = Bio::ToolBox->load_file($multi_file) or 
	die "can't load $multi_file!\n";
$MultiData->name(1, 'Start0'); # because it's actually 0-based, and this will trigger it
$MultiData->iterate( sub {
	my $row = shift;
	my $length = $row->length;
	my $key = $row->value(4); # this should the list column
	$name2space{$key} += $length;
	$total += $length;
});
$MultiData->save; # re-save with the new name.
undef $MultiData;

# write a spatial Venn file
my $VennData = Bio::ToolBox->new_data('File', 'BasePairs', 'Fraction');
foreach my $key (sort {$a cmp $b} keys %name2space) {
	$VennData->add_row( [$key, $name2space{$key}, 
		sprintf("%.3f", $name2space{$key} / $total) ]);
}
$VennData->add_comment("Coverage in bp for each category and fraction of total");
$VennData->save("$outfile\.spatialVenn.txt");



