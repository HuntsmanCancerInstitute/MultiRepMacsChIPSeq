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

use warnings;
use strict;
use English qw(-no_match_vars);
use Getopt::Long;
use File::Spec;
use File::Which;
use Bio::ToolBox 1.69;
use Bio::ToolBox::utility qw(format_with_commas);
use List::Util qw(min max uniqstr first);
use Statistics::Descriptive;

our $VERSION = 6.0;

# user variables
my $tool = which('bedtools');
my $outfile;
my $peak_basename;
my $min_peaks = 1;
my $genome_file;
my $merge_gap = 1;
my @files;
my $help;

# global variables
my $out_path;          # output directory
my %chromosomes;       # hash of wanted chromosomes
my @names;             # peak file base names
my @sortfiles;         # sorted peak file names
my %name2count;        # hash for total count of each peak file for later use
my $MergeData;         # merged peak Data object
my $total_bp   = 0;    # total basepairs over union
my $int_count  = 0;    # total count of intersections
my $peak_count = 0;    # number of accepted merged peaks
my %key2space;         # hash of basepair space for each intersection combination
my %key2count;         # hash of count for each intersection combination
my %current;           # reusable current intersection data hash

### Documentation
my $docs = <<DOC;

A script to intersect two or more peak interval files and generate a single 
union interval file, as well as numerous statistics describing the peaks and 
intersections. It uses the bedtools application for intersection calculations.

It will calculate a number of descriptive statistics for the input files, 
including count, sum, standard deviation, and quartiles of peak lengths.

It will run a multi-way intersection with the bedtools application and parse
the output into a merged interval bed file and an overlap summary statistics
file. A minimum number of overlapping peaks may be explicitly set to restrict
the merged peaks to those represented by multiple input peak files (default 1).

It will report the pairwise number of intersections and spatial overlap 
(Jaccard statistic) between all pairwise combinations of peak intervals. 

Results may be plotted using plot_peak_figures.R using the out basename as input.

Seven files will be written:
    output.bed                    the merged peaks in bed format
    output.matrix.txt             boolean intersection matrix for each merged peak
    output.jaccard.txt            pairwise Jaccard statistic (bp overlap) table
    output.n_intersection.txt     pairwise intersection count table
    output.multi.txt.gz           data file from bedtools multi-intersect tool 
    output.intersection.txt       intersection statistics for each peak combination
    output.lengthStats.txt        interval length statistics for each peak input

VERSION: $VERSION

USAGE: intersect_peaks.pl --out <basename> [options] <peak1> <peak2> ...

OPTIONS:
  -o --out <basename>      Provide the output basename
  -n --name <text>         Merged peak name (default output basename)
  -m --min <integer>       Minimum number of peak overlaps to merge (default 1)
  -a --gap <integer>       Maximum gap before merging neighboring peaks ($merge_gap bp)
  -g --genome <path>       Provide a genome file for sort consistency
  -b --bed <path>          Path to bedtools ($tool)
  -h --help                Print documentation
DOC

# Options
unless (@ARGV) {
	print $docs;
	exit;
}
GetOptions(
	'o|out=s'    => \$outfile,
	'n|name=s'   => \$peak_basename,
	'm|min=i'    => \$min_peaks,
	'a|gap=i'    => \$merge_gap,
	'g|genome=s' => \$genome_file,
	'b|bed=s'    => \$tool,
	'h|help!'    => \$help,
) or die "unrecognized option!\n";

### Main functions
check_options();
check_genome_file();
process_input_peak_files();
intersect_peaks();
write_venn_file();
calculate_jaccard();
cleanup();
exit;

################ Subroutines #######################################################

sub check_options {
	if ($help) {
		print $docs;
		exit;
	}
	unless ($outfile) {
		print "must provide output base name!\n";
		exit 1;
	}
	unless ($tool) {
		print "bedtools not in your PATH!\n";
		exit 1;
	}

	# strip any existing extension if provided
	$outfile =~ s/\. (?:bed|narrowPeak) (?:\.gz)? $//xi;
	my $out_base;
	( undef, $out_path, $out_base ) = File::Spec->splitpath($outfile);
	unless ($peak_basename) {
		$peak_basename = $out_base;
	}
	$out_path ||= '.';    # current directory

	# inputs
	push @files, @ARGV;
	unless ( scalar(@files) > 1 ) {
		print "must provide 2 or more input files!\n";
		exit 1;
	}
	if ( $genome_file and ( not -e $genome_file or not -r _ ) ) {
		print " Genome file '$genome_file' not present or readable!\n";
		exit 1;
	}
}

sub check_genome_file {

   # yeah, this may seem counter intuitive, but by sorting EVERY file in the SAME manner
   # I can avoid a bucket load of head aches because bedtools is a stickler for sort order
   # the problem is that I can't guarantee that the peaks files match the chromosome file
	return unless ($genome_file);

	my $Data = Bio::ToolBox->load_file(
		file     => $genome_file,
		noheader => 1
	) or die "can't read genome file '$genome_file'!";
	$Data->name( 0, 'chromosome' );
	$Data->name( 1, 'start' );
	$Data->gsort_data;
	$genome_file = File::Spec->catfile(
		$out_path,
		$Data->basename . ".sort.$PID." . $Data->extension
	);

	# we can't use standard file writing mechanism here, because it'll write headers
	# do it manually – this should be fixed in Bio::ToolBox v1.69
	# but leave in manual bit for now
	my $fh = Bio::ToolBox->write_file($genome_file)
		or die "unable to write '$genome_file'!";
	$Data->iterate(
		sub {
			my $row = shift;
			$fh->printf( "%s\t%s\n", $row->value(0), $row->value(1) );
			$chromosomes{ $row->value(0) } = $row->value(1);
		}
	);
	$fh->close;

}

sub process_input_peak_files {
	print " Collecting descriptive statistics on peak lengths...\n";

	# Peak length statistics file
	my $LengthStats = Bio::ToolBox->new_data(
		'File',        'Count',         'Sum',          'Mean',
		'StandardDev', 'Min',           'Percentile5',  'FirstQuartile',
		'Median',      'ThirdQuartile', 'Percentile95', 'Max'
	);
	$LengthStats->add_comment('Interval length statistics');

	# iterate through all the input files
	for my $file (@files) {

		# open file
		my $Data = Bio::ToolBox->load_file( file => $file );
		unless ($Data->number_rows > 0) {
			print " Unable to process file '$file', skipping\n";
			next;
		}
		unless ( $Data->bed ) {
			print "file '$file' does not appear to be BED-family format! skipping\n";
			next;
		}

		# process name
		my $name = $Data->basename;
		$name =~ s/\.rep [_\-\.] me (?: an | rge)$//x;   # remove rep_mean rep_merge
		$name =~ s/_peaks$//;
		if ( exists $name2count{$name} ) {
			print " basename '$name' provided more than once! skipping duplicates\n";
			next;
		}
		push @names, $name;
		printf "   %s has %d intervals\n", $name, $Data->last_row;

		# collect length numbers
		my $Stat = Statistics::Descriptive::Full->new();
		my @unwanted;
		$Data->iterate(
			sub {
				my $row = shift;
				if ( %chromosomes and not exists $chromosomes{ $row->seq_id } ) {
					push @unwanted, $row->row_index;
					next;
				}
				$Stat->add_data( $row->length );
			}
		);
		if (@unwanted) {
			printf "    excluding %d intervals not in genome file\n",
				scalar(@unwanted);
			$Data->delete_row(@unwanted);
		}

		# record stats
		$LengthStats->add_row(
			[
				$name,
				$Stat->count,
				$Stat->sum,
				sprintf( "%.0f", $Stat->mean ),
				sprintf( "%.0f", $Stat->standard_deviation ),
				$Stat->min,
				sprintf( "%.0f", $Stat->percentile(5) || $Stat->min ),
				sprintf( "%.0f", $Stat->quantile(1) ),
				sprintf( "%.0f", $Stat->quantile(2) ),                    # median
				sprintf( "%.0f", $Stat->quantile(3) ),
				sprintf( "%.0f", $Stat->percentile(95) || $Stat->max ),
				$Stat->max
			]
		);
		$name2count{$name} = $Stat->count;

		# sort the file just for sanity purposes
		$Data->gsort_data;
		my $sortfile = File::Spec->catfile(
			$out_path,
			$Data->basename . ".sort.$PID" . $Data->extension
		);
		$Data->save($sortfile);
		push @sortfiles, $sortfile;
	}

	# finish
	unless ( @sortfiles >= 2 ) {
		print " Need at least 2 valid files to proceed!\n";
		cleanup();
		exit 1;
	}
	$LengthStats->save("$outfile\.lengthStats.txt");
}

sub intersect_peaks {
	print " Calculating multi-way intersection....\n";

	my $multi_file = $outfile . '.multi.txt';
	my $command    = sprintf(
		"%s multiinter -header -names %s -i %s > %s",
		$tool,
		join( q( ), @names ),
		join( q( ), @sortfiles ), $multi_file
	);
	if ( system($command) ) {
		die "something went wrong! command:\n $command\n";
	}

	## parse intersections
	print " Parsing intersections....\n";
	my $MultiData = Bio::ToolBox->load_file($multi_file)
		or die "can't load $multi_file!\n";
	$MultiData->name( 1, 'Start0' );    # because it's actually 0-based
	$MultiData->add_comment(
		sprintf(
			"Output from bedtools multi-intersect tool between %s",
			join( ', ', @names )
		)
	);
	$MultiData->gsort_data;   # may not actually be sorted correctly!!!!!!!!!! 

	# initialize current merged peak
	%current = (
		chromo => $MultiData->value( 1, 0 ),
		start1 => 1,
		end1   => 1,
		names1 => q(),
		start2 => 1,
		end2   => 1,
		names2 => q(),
		score  => 0,
	);

	# output Data object for merged peaks
	$MergeData = Bio::ToolBox->new_data(
		'Chromosome', 'Start', 'End', 'Name', 'Score',
		@names
	);

	# iterate through the multi-intersections
	$MultiData->iterate( \&multi_intersect_callback );

	# save the last peak
	if ( $current{start1} != 1 and $current{end1} != 1 ) {
		process_merged_peak();
	}

	# re-save as compressed file
	unlink $multi_file;
	$multi_file .= '.gz';
	$MultiData->save($multi_file);

	## Write merged peak data file
	my $merge_peak_bed = $outfile . '.bed';
	my $merged_fh      = Bio::ToolBox->write_file($merge_peak_bed)
		or die "unable to write '$merge_peak_bed'!";
	$merged_fh->printf( "# merged peaks from %s\n", join( ', ', @names ) );
	$merged_fh->printf("# score represents sum of overlapping interval lengths\n");
	$MergeData->iterate(
		sub {
			my $row = shift;
			$merged_fh->printf(
				"%s\n",
				$row->bed_string(
					bed   => 5,
					score => $row->value(4),
				)
			);
		}
	);
	$merged_fh->close;
	my $merge_peak_file = $outfile . '.matrix.txt';
	$MergeData->delete_column( 0, 1, 2, 4 )
		;    # delete coordinates and count, no longer needed
	$MergeData->add_comment('Boolean matrix indicating overlapping source intervals');
	$MergeData->save($merge_peak_file);

	# summary
	printf "  %s independent intersection clusters were identified\n",
		format_with_commas($int_count);
	printf "  %s intervals were retained with %d or more overlaps\n",
		format_with_commas($peak_count), $min_peaks;
}

sub multi_intersect_callback {
	my $row = shift;

	# skip unwanted chromosomes
	return if ( %chromosomes and not exists $chromosomes{ $row->seq_id } );

	# process length
	my $length = $row->length;
	my $key    = join ',', sort { $a cmp $b } split( /,/, $row->value(4) );
	$total_bp += $length;
	$key2space{$key} += $length;

	# process the peak
	my $n = $row->value(3);    # reported number of overlapping peaks
	if ( $row->seq_id ne $current{chromo}
		or ( $row->start - $current{end1} ) > $merge_gap )
	{
		# we have moved to the next merged block
		# store the current valid interval
		if ( $current{start1} != 1 and $current{end1} != 1 ) {
			process_merged_peak();
		}

		# start the new one
		$current{chromo} = $row->seq_id;
		$current{start1} = $row->start;
		$current{end1}   = $row->end;
		$current{names1} = $key;
		if ( $n >= $min_peaks ) {
			$current{start2} = $row->start;
			$current{end2}   = $row->end;
			$current{score}  = ( $length * $n );
			$current{names2} = $key;
		}
		else {
			$current{start2} = 1;
			$current{end2}   = 1;
			$current{score}  = 0;
			$current{names2} = q();
		}
	}
	else {
		# same block, update end
		$current{end1} = $row->end;
		$current{names1} .= sprintf( ";%s", $key );
		if ( $n >= $min_peaks ) {
			if ( $current{start2} == 1 and $current{end2} == 1 ) {

				# start new block to keep
				$current{start2} = $row->start;
				$current{end2}   = $row->end;
				$current{score}  = ( $length * $n );
				$current{names2} = $key;
			}
			elsif ( $row->start - $current{end2} <= $merge_gap ) {

				# extend current block to keep
				$current{end2} = $row->end;
				$current{score} += ( $length * $n );
				$current{names2} .= sprintf( ";%s", $key );
			}
			else {
				# multiple sub peaks, record current one, then start new
				record_merged_peak();
				$current{start2} = $row->start;
				$current{end2}   = $row->end;
				$current{score}  = ( $length * $n );
				$current{names2} = $key;
			}
		}
	}
}

sub process_merged_peak {

	# this should work on global variables, nothing is passed

	# generate a representative key of all overlapping peak names
	my @bits       = sort { $a cmp $b } uniqstr( split /,|;/, $current{names1} );
	my $peakstring = join ',', @bits;

	# record intersection
	$key2count{$peakstring} += 1;
	$int_count++;

	# record merged interval if passes threshold
	if ( $current{start2} != 1 and $current{end2} != 1 ) {
		record_merged_peak();
	}
}

sub record_merged_peak {
	$peak_count++;
	my @bits = sort { $a cmp $b } uniqstr( split /,|;/, $current{names2} );
	my @data = (
		$current{chromo},
		$current{start2},
		$current{end2},
		sprintf( "%s.%d", $peak_basename, $peak_count ),
		$current{score},
	);

	# calculate boolean value whether each peak file overlaps this interval
	foreach my $n (@names) {
		if ( first { $_ eq $n } @bits ) {
			push @data, 1;
		}
		else {
			push @data, 0;
		}
	}
	$MergeData->add_row( \@data );
}

sub write_venn_file {
	print " Calculating intersection statistics....\n";
	my $VennData = Bio::ToolBox->new_data(
		'Peaks', 'BasePairs', 'BasePairFraction',
		'Count', 'CountFraction'
	);
	foreach my $key ( sort { $a cmp $b } keys %key2space ) {
		my @data =
			( $key, $key2space{$key}, sprintf( "%.3f", $key2space{$key} / $total_bp ) );
		if ( exists $key2count{$key} ) {
			push @data,
				( $key2count{$key}, sprintf( "%.3f", $key2count{$key} / $int_count ) );
		}
		else {
			push @data, ( 0, '0.000' );
		}
		$VennData->add_row( \@data );
	}
	$VennData->add_comment('Intersection statistics for each overlapping combination');
	$VennData->add_comment( sprintf( "Source intervals: %s", join( ',', @names ) ) );
	$VennData->add_comment("BasePairFraction is fraction of union length");
	$VennData->add_comment(
		"Count is number of independent clusters represented by the combination");
	my $venn_file = $outfile . '.intersection.txt';
	$VennData->save($venn_file);
}

sub calculate_jaccard {
	print " Calculating pairwise Jaccard overlap....\n";

	# output for intersection data
	my $JaccardData = Bio::ToolBox->new_data( 'File', @names );
	foreach my $n (@names) {
		$JaccardData->add_row( [ $n, ( map {0} @names ) ] );
	}
	$JaccardData->add_comment(
		"Pairwise fraction overlap of the union between each source file");
	my $IntersectionData = Bio::ToolBox->new_data( 'File', @names );
	foreach my $n (@names) {
		$IntersectionData->add_row( [ $n, ( map {0} @names ) ] );
	}
	$IntersectionData->add_comment(
		"Number of pairwise intersections between each source file");
	my %name2i;
	{
		my $i = 1;
		%name2i = map { $_ => $i++ } @names;
	}

	# bedtools jaccard base command
	my $jaccard_cmdbase = "$tool jaccard ";
	$jaccard_cmdbase .= "-g $genome_file " if $genome_file;
	my $jaccard_warn;

	# loop through files
	my %reciprocal;
	foreach my $f1 ( 0 .. $#names ) {

		foreach my $f2 ( 0 .. $#names ) {

			# skip those we don't need to do
			if ( exists $reciprocal{"$f1\_$f2"} ) {

				# don't need to do reciprocal
				next;
			}
			elsif ( $f1 == $f2 ) {

				# self intersection - we know the answer to this
				my $i = $name2i{ $names[$f1] };
				$JaccardData->value( $i, $i, '1.0000' );
				$IntersectionData->value( $i, $i, $name2count{ $names[$f1] } );
				$reciprocal{"$f2\_$f1"} = 1;
				next;
			}

			# coordinates
			my $x = $name2i{ $names[$f1] };
			my $y = $name2i{ $names[$f2] };

			my $result = qx($jaccard_cmdbase -a $sortfiles[$f1] -b $sortfiles[$f2] 2>&1);
			if ( $result =~ /error/i ) {

				# there was some sort of error - the most likely is lexicographic
				# the user needs a genome file because the bed files are either not
				# sorted properly or intervals are not present on all the chromosomes
				printf "   jaccard between %s and %s failed!!\n", $names[$f1],
					$names[$f2];
				$jaccard_warn = $result;
				$JaccardData->value( $x, $y, '.' );
				$JaccardData->value( $y, $x, '.' );
				$IntersectionData->value( $x, $y, '.' );
				$IntersectionData->value( $y, $x, '.' );
			}
			else {
				# likely worked
				printf "   jaccard between %s and %s succeeded\n", $names[$f1],
					$names[$f2];
				my @lines = split( /\n/, $result );
				my ( undef, undef, $jaccard, $number ) = split( /\s+/, $lines[1] );
				$JaccardData->value( $x, $y, sprintf( "%.4f", $jaccard ) );
				$JaccardData->value( $y, $x, sprintf( "%.4f", $jaccard ) );
				$IntersectionData->value( $x, $y, $number );
				$IntersectionData->value( $y, $x, $number );
			}

			# done
			$reciprocal{"$f2\_$f1"} = 1;
		}
	}
	if ($jaccard_warn) {
		print
"  Jaccard calculations are incomplete!!! This was the last error:\n$jaccard_warn";
	}

	# write files
	$JaccardData->save("$outfile\.jaccard.txt");
	$IntersectionData->save("$outfile\.n_intersection.txt");
}

sub cleanup {
	unlink @sortfiles;
	if ($genome_file) {
		unlink $genome_file;
	}
}

