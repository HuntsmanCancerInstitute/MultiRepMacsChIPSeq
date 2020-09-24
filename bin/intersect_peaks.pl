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

my $VERSION = 3;

# variables
my $tool = which('bedtools');
my $outfile;
my $peak_basename;
my $genome_file;
my @files;
my $help;


### Documentation
my $docs = <<DOC;

A script to intersect two or more peak files. This is a wrapper around the bedtools 
program.

It will first merge all of the peak files into a single representative bed file.
Peak intervals will be renamed to the given name. 

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
    --name text             Provide text to rename the merged peaks
    --genome path           Provide a genome file for consistency
    --bed path              Path to bedtools ($tool)
    --help                  Print documentation
DOC




### Options
unless (@ARGV) {
	print $docs;
	exit;
}
GetOptions(
	'out=s'             => \$outfile,
	'name=s'            => \$peak_basename,
	'genome=s'          => \$genome_file,
	'bed=s'             => \$tool,
	'help!'             => \$help,
) or die "unrecognized option!\n";

# check options
if ($help) {
	print $docs;
	exit;
}
die "must provide output base name!\n" unless $outfile;
die "bedtools not in your PATH!\n" unless $tool;
$outfile =~ s/\.(?:bed|narrowPeak)$//i; # strip any existing extension if provided
unless ($peak_basename) {
	# use the outfile basename
	(undef, undef, $peak_basename) = File::Spec->splitpath($outfile);
}

# inputs
my @files = @ARGV;
die "must provide 2 or more files!\n" unless scalar(@files) > 1;
my @names;
foreach my $f (@files) {
	# strip path and extension
	unless (-e $f and -r _ ) {
		die " Input file '$f' not present or readable!\n";
	}
	my (undef, $dir, $name) = File::Spec->splitpath($f);
	$name =~ s/\.(?:bed|narrowpeak|gappedpeak|broadpeak)(?:\.gz)?$//i;
	push @names, $name;
}
if ($genome_file and (not -e $genome_file or not -r _ )) {
	die " Genome file '$genome_file' not present or readable!\n";
}



### Merge the peaks
print " Merging peak files....\n";
my $merge_peak_file = $outfile . '.bed';
my $command = sprintf("cat %s | %s sort -i - | %s merge -c 5 -o mean -i - > %s",
	join(" ", @files), $tool, $tool, $merge_peak_file);
if (system($command)) {
	die "something went wrong! command:\n $command\n";
}

# Clean up the merged peaks
if (-e $merge_peak_file and -s _ ) {
	my $Data = Bio::ToolBox->load_file($merge_peak_file) or 
		die "unable to open $merge_peak_file!";
	
	# sort the file correctly, not asciibetically
	$Data->gsort_data;
	
	# rename
	$Data->iterate( sub {
		my $row = shift;
		$row->value(3, sprintf("%s.%d", $peak_basename, $row->row_index));
	});
	$Data->save;
}



### Calculate Jaccard statistic
print " Calculating Jaccard overlap....\n";
my $JaccardData = Bio::ToolBox->new_data('File', @names);
my $IntersectionData = Bio::ToolBox->new_data('File', @names);
my $jaccard_cmdbase = "$tool jaccard ";
$jaccard_cmdbase .= "-g $genome_file " if -e $genome_file;
my $jaccard_warn;

for my $f1 (0..$#names) {
	my $j = $JaccardData->add_row();
	my $i = $IntersectionData->add_row();
	
	$JaccardData->value($j, 0, $names[$f1]);
	$IntersectionData->value($i, 0, $names[$f1]);
	
	for my $f2 (0..$#names) {
		my $result = qx($jaccard_cmdbase -a $files[$f1] -b $files[$f2] 2>&1);
		if ($result =~ /ERROR/) {
			# there was some sort of error - the most likely is lexicographic 
			# the user needs a genome file because the bed files are either not 
			# sorted properly or intervals are not present on all the chromosomes
			printf "   jaccard between %s and %s failed!!\n", $names[$f1], $names[$f2];
			$jaccard_warn = $result;
			$JaccardData->value($j, $f2 + 1, '.');
			$IntersectionData->value($i, $f2 + 1, '.');
		}
		else {
			# likely worked
			printf "   jaccard between %s and %s succeeded\n", $names[$f1], $names[$f2];
			my @lines = split("\n", $result);
			my (undef, undef, $jaccard, $number) = split(/\s+/, $lines[1]);
			$JaccardData->value($j, $f2 + 1, $jaccard);
			$IntersectionData->value($i, $f2 + 1, $number);
		}
	}
}
if ($jaccard_warn) {
	print "  Jaccard calculations are incomplete!!! This was the last error:\n$jaccard_warn";
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



