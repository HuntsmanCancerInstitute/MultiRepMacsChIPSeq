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
use List::Util qw(min max uniqstr);
use Statistics::Descriptive;

my $VERSION = 5.0;

# variables
my $tool = which('bedtools');
my $outfile;
my $peak_basename;
my $genome_file;
my $merge_gap = 1;
my @files;
my $help;


### Documentation
my $docs = <<DOC;

A script to intersect two or more peak interval files and generate a single 
union interval file, as well as numerous statistics describing the peaks and 
intersections. It uses the bedtools application for intersection calculations.

It will first calculate a number of descriptive statistics for the input files, 
including count, sum, standard deviation, and quartiles of peak lengths.

It will then run a multi-way intersection with the bedtools application, and 
the output parsed into a union interval bed file and a summary statistics.

It will also report the pairwise number of intersections and spatial overlap 
(Jaccard statistic) between all peak intervals. 


Seven files will be written:
    output.bed                    the merged peaks in bed format
    output.txt                    the merged peaks in text format
    output.jaccard.txt            pairwise Jaccard statistic (bp overlap) table
    output.n_intersection.txt     pairwise intersection count table
    output.multi.txt              data file from bedtools multi-intersect tool 
    output.intersection.txt       intersection statistics for each peak combination
    output.lengthStats.txt        interval length statistics for each peak input

VERSION: $VERSION

USAGE: intersect_peaks.pl --out <basename> peak1.narrowPeak peak2.narrowPeak ....

OPTIONS:
    --out basename          Provide the output basename
    --name text             Merged peak name (default out basename)
    --gap integer           Maximum gap before merging neighboring peaks ($merge_gap bp)
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
	'gap=i'             => \$merge_gap,
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



### Check genome file
# yeah, this may seem counter intuitive, but by sorting EVERY file in the SAME manner
# I can avoid a bucket load of head aches because bedtools is a stickler for sort order
# the problem is that I can't guarantee that the peaks files match the chromosome file
my %chromosomes;
if ($genome_file) {
	my $Data = Bio::ToolBox->load_file(
		file     => $genome_file,
		noheader => 1
	) or die "can't read genome file '$genome_file'! $!";
	$Data->name(0, 'chromosome');
	$Data->name(1, 'start');
	$Data->gsort_data;
	$genome_file = File::Spec->catfile($out_path, $Data->basename . ".sort.$$." . $Data->extension);
	# we can't use standard file writing mechanism here, because it'll write headers
	# do it manually – this should be fixed in Bio::ToolBox v1.69 but leave in manual bit for now
	my $fh = Bio::ToolBox->write_file($genome_file) or 
		die "unable to write $genome_file! $!";
	$Data->iterate( sub {
		my $row = shift;
		$fh->printf("%s\t%s\n", $row->value(0), $row->value(1));
		$chromosomes{ $row->value(0) } = $row->value(1);
	});
	$fh->close;
}


### Interval length statistics
print " Collecting descriptive statistics on peak lengths...\n";
my @names;
my $LengthStats = Bio::ToolBox->new_data('File', 'Count', 'Sum', 
	'Mean', 'StandardDev', 'Min', 'Percentile5','FirstQuartile', 'Median', 
	'ThirdQuartile', 'Percentile95', 'Max');
my @sortfiles;
my %name2count; # hash for total count of each peak file for later use
for my $file (@files) {
	
	# open file
	my $Data = Bio::ToolBox->load_file(file => $file)
		or die "unable to open '$file' $!"; 
	unless ($Data->bed) {
		die "file '$file' does not appear to be a BED file format!\n";
	}
	printf "   %s has %d intervals\n", $file, $Data->last_row;
	
	# remember stuff about this file
	push @names, $Data->basename;
	
	# collect length numbers
	my $Stat = Statistics::Descriptive::Full->new();
	my @unwanted;
	$Data->iterate( sub {
		my $row = shift;
		if (%chromosomes and not exists $chromosomes{ $row->seq_id }) {
			push @unwanted, $row->row_index;
			next;
		}
		$Stat->add_data( $row->length );
	} );
	if (@unwanted) {
		printf "    excluding %d intervals not in genome file for intersection\n", scalar(@unwanted);
		$Data->delete_row(@unwanted);
	}
	
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
	$name2count{ $Data->basename } = $Stat->count;
	
	# sort the file just for sanity purposes
	$Data->gsort_data;
	my $sortfile = File::Spec->catfile($out_path, $Data->basename . ".sort.$$" . $Data->extension);
	$Data->save($sortfile);
	push @sortfiles, $sortfile;
}
$LengthStats->save("$outfile\.lengthStats.txt");





### Multiple intersection
print " Calculating multi-way intersection....\n";
my $multi_file = $outfile . '.multi.txt';
my $command = sprintf("%s multiinter -header -names %s -i %s > %s", 
	$tool, join(' ', @names), join(' ', @sortfiles), $multi_file);
if (system($command)) {
	die "something went wrong! command:\n $command\n";
}

## parse intersections
print " Parsing intersections....\n";
my $MultiData = Bio::ToolBox->load_file($multi_file) or 
	die "can't load $multi_file!\n";
$MultiData->name(1, 'Start0'); # because it's actually 0-based, and this will trigger it
$MultiData->add_comment(
	sprintf("Output from bedtools multi-intersect tool between %s", 
	join(', ', @names))
);

# initialize counters
my $total_bp = 0;
my $total_count = 0;
my %key2space;
my %key2count;
my %current = (
	chromo => $MultiData->value(1,0),
	start  => 1,
	end    => 1,
	names  => '',
);

# output for merged peak
my $MergeData = Bio::ToolBox->new_data('Chromosome', 'Start', 'End', 'Name', 'NumberIntervals', 'Peaks');
$MergeData->add_comment(sprintf("Merged '%s' peaks from %s", $peak_basename, 
	join(', ', @names)));

# iterate through the multi-intersections
$MultiData->iterate( sub {
	my $row = shift;
	
	# skip unwanted chromosomes
	next if (%chromosomes and not exists $chromosomes{ $row->seq_id });
	
	# process length
	my $length = $row->length;
	my @keys = sort {$a cmp $b} split(',', $row->value(4));
	my $key = join(',', @keys);
	$key2space{$key} += $length;
	$total_bp += $length;
	
	# check current window
	if (
		$row->seq_id ne $current{chromo} or 
		($row->start - $current{end}) > $merge_gap
	) {
		# we have moved to the next merged block
		# store the current one
		unless ($current{start} == 1 and $current{end} == 1) {
			# but not if it's the fake startup data!
			process_merged_peak();
		}
		
		# start the new one
		$current{chromo} = $row->seq_id;
		$current{start}  = $row->start;
		$current{end}    = $row->end;
		$current{names}  = $key;
	}
	else {
		# same block, update end
		$current{end} = $row->end;
		if ($current{names}) {
			$current{names} .= sprintf(";%s", $key);
		}
		else {
			$current{names}  = $key;
		}
	}
} );

# save the last peak
process_merged_peak();

# re-save
$MultiData->save; 
undef $MultiData;




## Write an intersection Venn file
my $VennData = Bio::ToolBox->new_data('Peaks', 'BasePairs', 'BasePairFraction', 
	'Count', 'CountFraction');
foreach my $key (sort {$a cmp $b} keys %key2space) {
	my @data = ($key, $key2space{$key}, sprintf("%.3f", $key2space{$key} / $total_bp) );
	if (exists $key2count{$key} ) {
		push @data, ($key2count{$key}, sprintf("%.3f", $key2count{$key} / $total_count));
	}
	else {
		push @data, (0, '0.000');
	}
	$VennData->add_row(\@data);
}
$VennData->add_comment("Intersection coverage in bp and fraction of total for each input combination");
$VennData->add_comment("Intersection count and fraction of total for each input combination");
my $venn_file = $outfile . '.intersection.txt';
$VennData->save($venn_file);



## Write merged peak data file
my $merge_peak_bed = $outfile . '.bed';
my $merged_fh = Bio::ToolBox->write_file($merge_peak_bed) or 
	die "unable to write $merge_peak_bed! $!";
$MergeData->iterate( sub {
	my $row = shift;
	$merged_fh->printf("%s\n", $row->bed_string(bed => 4));
} );
$merged_fh->close;
my $merge_peak_file = $outfile . '.txt';
$MergeData->save($merge_peak_file);
printf "  wrote %d intersected peaks to $merge_peak_bed and $merge_peak_file\n", $total_count;
printf "  wrote intersection statistics to $venn_file\n";








### Calculate Jaccard statistic
print " Calculating Jaccard overlap....\n";

# output for intersection data
my $JaccardData = Bio::ToolBox->new_data('File', @names);
foreach my $n (@names) {
	$JaccardData->add_row([$n, (map {0} @names)]);
}
$JaccardData->add_comment("Pairwise fraction overlap of the union between each input file");
my $IntersectionData = Bio::ToolBox->new_data('File', @names);
foreach my $n (@names) {
	$IntersectionData->add_row([$n, (map {0} @names)]);
}
$IntersectionData->add_comment("Number of pairwise intersections between each input file");
my %name2i;
{
	my $i = 1;
	%name2i = map {$_ => $i++} @names;
}

# bedtools jaccard base command
my $jaccard_cmdbase = "$tool jaccard ";
$jaccard_cmdbase .= "-g $genome_file " if $genome_file;
my $jaccard_warn;

# loop through files
my %reciprocal;
foreach my $f1 (0..$#names)  {
	
	foreach my $f2 (0..$#names) {
		
		# skip those we don't need to do
		if (exists $reciprocal{"$f1\_$f2"}) {
			# don't need to do reciprocal
			next;
		}
		elsif ($f1 == $f2) {
			# self intersection - we know the answer to this
			my $i = $name2i{$names[$f1]};
			$JaccardData->value($i, $i, '1.0000');
			$IntersectionData->value($i, $i, $name2count{$names[$f1]});
			$reciprocal{"$f2\_$f1"} = 1; 
			next;
		}
		
		# coordinates
		my $x = $name2i{ $names[$f1] };
		my $y = $name2i{ $names[$f2] };
		
		my $result = qx($jaccard_cmdbase -a $sortfiles[$f1] -b $sortfiles[$f2] 2>&1);
		if ($result =~ /error/i) {
			# there was some sort of error - the most likely is lexicographic 
			# the user needs a genome file because the bed files are either not 
			# sorted properly or intervals are not present on all the chromosomes
			printf "   jaccard between %s and %s failed!!\n", $names[$f1], $names[$f2];
			$jaccard_warn = $result;
			$JaccardData->value($x,$y, '.');
			$JaccardData->value($y,$x, '.');
			$IntersectionData->value($x,$y, '.');
			$IntersectionData->value($y,$x, '.');
		}
		else {
			# likely worked
			printf "   jaccard between %s and %s succeeded\n", $names[$f1], $names[$f2];
			my @lines = split("\n", $result);
			my (undef, undef, $jaccard, $number) = split(/\s+/, $lines[1]);
			$JaccardData->value($x,$y, sprintf("%.4f", $jaccard));
			$JaccardData->value($y,$x, sprintf("%.4f", $jaccard));
			$IntersectionData->value($x,$y, $number);
			$IntersectionData->value($y,$x, $number);
		}
		
		# done
		$reciprocal{"$f2\_$f1"} = 1; 
	}
}
if ($jaccard_warn) {
	print "  Jaccard calculations are incomplete!!! This was the last error:\n$jaccard_warn";
}
$JaccardData->save("$outfile\.jaccard.txt");
$IntersectionData->save("$outfile\.n_intersection.txt");






### clean up
unlink @sortfiles;
if ($genome_file) {
	unlink $genome_file;
}



sub process_merged_peak {
	# this should work on global variables, nothing is passed
	
	# generate a representative key of all overlapping peak names
	my $peakstring = join(',',
		sort {$a cmp $b} 
		uniqstr( split(/,|;/, $current{names}) )
	);
	
	# add merged peak data
	$key2count{$peakstring} += 1;
	$total_count++;
	$MergeData->add_row( [
		$current{chromo},
		$current{start},
		$current{end},
		sprintf("%s.%d", $peak_basename, $total_count),
		scalar(split(';', $current{names})),
		$peakstring
	] );
	
}

