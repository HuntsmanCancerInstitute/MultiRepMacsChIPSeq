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
use Getopt::Long;
use List::Util qw(max);
use Bio::ToolBox 1.69;

our $VERSION = 1;

### Documentation
my $docs = <<DOC;

A script to fill in missing score values in a Peak file.

Running MACS2 'bdgpeakcall' or 'bdgbroadcall' functions will write an
appropriate Peak output file, but will have missing score values,
specifically: 

    signalValue  - Fold Enrichment at summit or across broad region
    pValue       - P-Value at summit or across broad region
    qValue       - Q-Value at summit or across broad region 

One or more of these values may be updated by providing the appropriate 
score track file (currently only bigWig tracks are supported). 

Additionally, the file is automatically re-sorted with a sane chromosome 
sort, and names may be adjusted if desired, using the provided base text
that will be appended with an increasing number. Scores can also optionally 
be scaled to the UCSC-recommended range of 0-1000. 

For narrowPeak files, a summit bed file may be exported, if desired.

USAGE: update_peak_file.pl -i <peakfile> [options]

OPTIONS:
   -i --in <file>         Input narrowPeak file
   -o --out <file>        Output file, default input
   -n --name <text>       Base text to optionally rename the peaks
   -e --enrich <file>     The enrichment score bigWig file, e.g. log2FE
   -p --pvalue <file>     The P-Calue bigWig score file
   -q --qvalue <file>     The Q-Value bigWig score file
   -r --norm              Normalize the Score column to range 0-1000
   -u --summit            Export a summit Bed file from narrowPeak
   -h --help              Print documentation

DOC

# variables
my $infile;
my $outfile;
my $name_base;
my $enrich_file;
my $pval_file;
my $qval_file;
my $normalize;
my $export_summit;
my $help;
my $base;

unless (@ARGV) {
	print $docs;
	exit;
}

### Options
GetOptions(
	'i|input=s'  => \$infile,
	'o|out=s'    => \$outfile,
	'n|name=s'   => \$name_base,
	'e|enrich=s' => \$enrich_file,
	'p|pvalue=s' => \$pval_file,
	'q|qvalue=s' => \$qval_file,
	'r|norm!'    => \$normalize,
	'u|summit!'  => \$export_summit,
	'help!'      => \$help,
) or die "unrecognized option! See help\n";

# Check options
if ($help) {
	print $docs;
	exit;
}
unless ($infile) {
	print " No input narrowPeak file provided!\n";
	exit 1;
}
if ( $enrich_file and $enrich_file =~ m/log (\d{1,2})/xi ) {
	$base = $1;
}

# Open input
my $Data = Bio::ToolBox->load_file($infile);
unless ($Data) {
	print " Unable to load '$infile'!\n";
	exit 1;
}
unless ( $Data->bed and $Data->format =~ /peak/i ) {
	print " Input file is not a Peak format!\n";
	exit 1;
}

# Process
$Data->gsort_data;
my $max_score = collect_max_score();
print_summary();
my $Summit;
if ($export_summit) {
	$Summit = Bio::ToolBox->new_data(qw(Chromosome Start0 End Name));
	$Summit->bed(4);
}
if ( $Data->format eq 'narrowPeak' ) {
	$Data->iterate( \&narrowpeak_callback );
}
elsif ( $Data->format eq 'gappedPeak' ) {
	$Data->iterate( \&gappedpeak_callback );
}
elsif ( $Data->format eq 'broadPeak' ) {
	$Data->iterate( \&broadpeak_callback );
}
else {
	die sprintf( " unrecognized file format '%s'!", $Data->format );
}

# Finish
if ($outfile) {
	$Data->add_file_metadata($outfile);
}
my $s = $Data->save;
if ($s) {
	print " Updated values and wrote '$s'\n";
}
else {
	print " Problems saving '$s'!!!\n";
}
if ($export_summit) {
	my $sfile = $Data->path . $Data->basename . '.summit.bed';
	$Summit->add_file_metadata($sfile);
	my $sf = $Summit->save;
	if ($sf) {
		print " Saved summit file as '$sf'\n";
	}
	else {
		print " Problems saving '$sf'!!!\n";
	}
}

########### Subroutines ############
sub collect_max_score {
	my @scores = $Data->column_values(4);
	shift @scores;    # discard name
	return max(@scores);
}

sub print_summary {
	printf " Loaded '%s' with %d intervals\n", $infile, $Data->last_row;
	if ($name_base) {
		printf "  - Naming intervals with '%s'\n", $name_base;
	}
	if ($max_score) {
		printf "  - Scaling Score based on maximum value %s\n", $max_score;
	}
	if ($enrich_file) {
		printf "  - Adding SignalValue from file '%s'\n", $enrich_file;
		if ($base) {
			printf "    De-logging with base %d\n", $base;
		}
	}
	if ($pval_file) {
		printf "  - Adding P-Value from file '%s'\n", $pval_file;
	}
	if ($qval_file) {
		printf "  - Adding Q-Value from file '%s'\n", $qval_file;
	}
}

sub narrowpeak_callback {
	my $row = shift;

	# check name
	if ($name_base) {
		my $n = sprintf "%s.%d", $name_base, $row->row_index;
		$row->value( 3, $n );
	}

	# peak coordinate
	my $p = $row->peak;
	unless ($p) {
		die sprintf " No peak coordinate for row %d!", $row->row_index;
	}

	# check enrichment signal value
	if ($enrich_file) {
		my $v = $row->get_score(
			dataset => $enrich_file,
			start   => $p,
			stop    => $p
		);
		if ( defined $base ) {
			$v = $base**$v;
		}
		$row->value( 6, sprintf( "%.3f", $v ) );
	}

	# check p-value
	if ($pval_file) {
		my $v = $row->get_score(
			dataset => $pval_file,
			start   => $p,
			stop    => $p
		);
		$row->value( 7, sprintf( "%.3f", $v ) );
	}

	# check q-value
	if ($qval_file) {
		my $v = $row->get_score(
			dataset => $qval_file,
			start   => $p,
			stop    => $p
		);
		$row->value( 8, sprintf( "%.3f", $v ) );
	}

	# normalize
	if ($normalize) {
		my $v = $row->value(4);
		$row->value( 4, sprintf( "%.0f", ( $v / $max_score ) * 1000 ) );
	}

	# export summit
	if ($export_summit) {
		$Summit->add_row(
			[
				$row->seq_id,
				$p - 1,
				$p,
				$row->value(3)
			]
		);
	}
}

sub gappedpeak_callback {
	my $row = shift;

	# check name
	if ($name_base) {
		my $n = sprintf "%s.%d", $name_base, $row->row_index;
		$row->value( 3, $n );
	}

	# check enrichment signal value
	if ($enrich_file) {
		my $v = $row->get_score(
			dataset => $enrich_file,
			method  => 'mean',
		);
		if ( defined $base ) {
			$v = $base**$v;
		}
		$row->value( 12, sprintf( "%.3f", $v ) );
	}

	# check p-value
	if ($pval_file) {
		my $v = $row->get_score(
			dataset => $pval_file,
			method  => 'mean',
		);
		$row->value( 13, sprintf( "%.3f", $v ) );
	}

	# check q-value
	if ($qval_file) {
		my $v = $row->get_score(
			dataset => $qval_file,
			method  => 'mean',
		);
		$row->value( 14, sprintf( "%.3f", $v ) );
	}

	# normalize
	if ($normalize) {
		my $v = $row->value(4);
		$row->value( 4, sprintf( "%.0f", ( $v / $max_score ) * 1000 ) );
	}
}

sub broadpeak_callback {
	my $row = shift;

	# check name
	if ($name_base) {
		my $n = sprintf "%s.%d", $name_base, $row->row_index;
		$row->value( 3, $n );
	}

	# check enrichment signal value
	if ($enrich_file) {
		my $v = $row->get_score(
			dataset => $enrich_file,
			method  => 'mean',
		);
		if ( defined $base ) {
			$v = $base**$v;
		}
		$row->value( 6, sprintf( "%.3f", $v ) );
	}

	# check p-value
	if ($pval_file) {
		my $v = $row->get_score(
			dataset => $pval_file,
			method  => 'mean',
		);
		$row->value( 7, sprintf( "%.3f", $v ) );
	}

	# check q-value
	if ($qval_file) {
		my $v = $row->get_score(
			dataset => $qval_file,
			method  => 'mean',
		);
		$row->value( 8, sprintf( "%.3f", $v ) );
	}

	# normalize
	if ($normalize) {
		my $v = $row->value(4);
		$row->value( 4, sprintf( "%.0f", ( $v / $max_score ) * 1000 ) );
	}
}

