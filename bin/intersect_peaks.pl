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
use List::Util qw(min max);
use Statistics::Descriptive;

my $VERSION = 4.1;

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

It will first calculate a number of descriptive statistics for the interval lengths 
for each input file, including count, sum, standard deviation, and quartiles.

It will merge all of the peak files into a single representative bed file.
Peak intervals will be renamed to the given name. 

It will run the bedtools Jaccard statistic pairwise across all of the peak 
files and write out a merged table of the results. The Jaccard statistic measures 
the amount of spatial overlap between two peak files (intersection/union) reported 
as a fraction between 0 and 1.

Finally, it will run bedtools multiinter tool to perform a multi-way intersection 
and calculate the intervals for each category of overlap. This is parsed into a 
summary file suitable for drawing a Venn diagram based on spatial overlap.

Six files will be written:
    basename.bed                    the merged peaks
    basename.jaccard.txt            the Jaccard results in a table
    basename.n_intersection.txt     the number of intersections in a table
    basename.multi.txt              data file from multi-intersection 
    basename.spatialVenn.txt        summary of spatial overlap for each category
    basename.lengthStats.txt        interval length statistics for each file

VERSION: $VERSION

USAGE: intersect_peaks.pl --out <basename> peak1.narrowPeak peak2.narrowPeak ....

OPTIONS:
    --out basename          Provide the output basename
    --name text             Provide text to rename the merged peaks
    --genome path           Provide a genome file for sort consistency
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
my (undef, $out_path, $out_base) = File::Spec->splitpath($outfile);
unless ($peak_basename) {
	# use the outfile basename
	$peak_basename = $out_base;
}
$out_path ||= '.'; # current directory

# inputs
my @files = @ARGV;
die "must provide 2 or more files!\n" unless scalar(@files) > 1;
if ($genome_file and (not -e $genome_file or not -r _ )) {
	die " Genome file '$genome_file' not present or readable!\n";
}


### Interval length statistics
print " Collecting descriptive statistics on peak lengths...\n";
my @numcols; # remember how many columns each input file has
my @names;
my $LengthStats = Bio::ToolBox->new_data('File', 'Count', 'Sum', 
	'Mean', 'StandardDev', 'Min', 'Percentile5','FirstQuartile', 'Median', 
	'ThirdQuartile', 'Percentile95', 'Max');
my @sortfiles;
for my $file (@files) {
	
	# open file
	print "   $file\n";
	my $Data = Bio::ToolBox->load_file(file => $file)
		or die "unable to open '$file' $!"; 
	unless ($Data->bed) {
		die "file '$file' does not appear to be a BED file format!\n";
	}
	
	# remember stuff about this file
	push @names, $Data->basename;
	push @numcols, $Data->number_columns;
	
	# collect length numbers
	my $Stat = Statistics::Descriptive::Full->new();
	$Data->iterate( sub {
		my $feature = shift;
		$Stat->add_data( $feature->length );
	} );
	
	# record stats
	$LengthStats->add_row( [
		$Data->basename,
		$Stat->count,
		$Stat->sum,
		sprintf("%.0f", $Stat->mean),
		sprintf("%.0f", $Stat->standard_deviation),
		$Stat->min,
		sprintf("%.0f", $Stat->percentile(5)),
		sprintf("%.0f", $Stat->quantile(1)),
		sprintf("%.0f", $Stat->quantile(2)), # median
		sprintf("%.0f", $Stat->quantile(3)),
		sprintf("%.0f", $Stat->percentile(95)),
		$Stat->max
	] );
	
	# sort the file just for sanity purposes
	$Data->gsort_data;
	my $sortfile = File::Spec->catfile($out_path, $Data->basename . '.sort' . $Data->extension);
	$Data->save($sortfile);
	push @sortfiles, $sortfile;
}
$LengthStats->save("$outfile\.lengthStats.txt");


### Sort genome file
# yeah, this may seem counter intuitive, but by sorting EVERY file in the SAME manner
# I can avoid a bucket load of head aches because bedtools is a stickler for sort order
# the problem is that I can't guarantee that the peaks files match the chromosome file
if ($genome_file) {
	my $Data = Bio::ToolBox->load_file(
		file     => $genome_file,
		noheader => 1
	) or die "can't read genome file '$genome_file'! $!";
	$Data->name(0, 'chromosome');
	$Data->name(1, 'start');
	$Data->gsort_data;
	$genome_file = File::Spec->catfile($out_path, $Data->basename . '.sort' . $Data->extension);
	# we can't use standard file writing mechanism here, because it'll write headers
	# do it manually
	my $fh = Bio::ToolBox->write_file($genome_file) or 
		die "unable to write $genome_file! $!";
	$Data->iterate( sub {
		my $row = shift;
		$fh->printf("%s\t%s\n", $row->value(0), $row->value(1));
	});
	$fh->close;
}


### Merge the peaks
print " Merging peak files....\n";
my $merge_peak_file = $outfile . '.bed';
my $command = sprintf("cat %s | ", join(" ", @sortfiles));
if (min(@numcols) != max(@numcols)) {
	# unequal columns, must trim excess columns
	$command .= sprintf("cut -f%d-%d | ", 1, min(@numcols));
}
$command .= sprintf("%s sort ", $tool);
if ($genome_file) {
	$command .= "-faidx $genome_file ";
}
$command .= sprintf("-i - | %s merge ", $tool);
if (min(@numcols) >= 5) {
	# record an average score
	$command .= "-c 5 -o mean ";
}
$command .= sprintf("-i - > %s", $merge_peak_file);
# print "   Executing $command\n";
if (system($command)) {
	die "something went wrong! command:\n $command\n";
}

# Clean up the merged peaks
if (-e $merge_peak_file and -s _ ) {
	my $Data = Bio::ToolBox->load_file($merge_peak_file) or 
		die "unable to open $merge_peak_file!";
	
	# sort the file correctly, not asciibetically - necessary????
	$Data->gsort_data;
	
	# rename intervals
	if ($Data->number_columns == 3) {
		$Data->add_column('Name'); 
	}
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
$jaccard_cmdbase .= "-g $genome_file " if $genome_file;
my $jaccard_warn;

for my $f1 (0..$#names) {
	my $j = $JaccardData->add_row();
	my $i = $IntersectionData->add_row();
	
	$JaccardData->value($j, 0, $names[$f1]);
	$IntersectionData->value($i, 0, $names[$f1]);
	
	for my $f2 (0..$#names) {
		my $result = qx($jaccard_cmdbase -a $sortfiles[$f1] -b $sortfiles[$f2] 2>&1);
		if ($result =~ /error/i) {
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
	$tool, join(' ', @names), join(' ', @sortfiles), $multi_file);
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



### clean up
unlink @sortfiles;
if ($genome_file) {
	unlink $genome_file;
}

