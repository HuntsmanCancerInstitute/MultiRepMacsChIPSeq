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
use English qw(-no_match_vars);
use File::Spec;
use Getopt::Long;
use Bio::ToolBox 1.65;

our $VERSION = 2;

### Documentation
my $docs = <<DOC;

A script to convert Macs2 narrowPeak and gappedPeak files into simpler bed files. 
Specifically, discard unused p-value and q-value columns after running Macs2 
bdgpeakcall or bdgbroadcall. 

This script will also resort the peaks properly (numerically, as opposed to 
alphabetically) and rename the peaks.

For narrowPeak files, one or two files may be written:
	a 5-column BED file of peak intervals, including original score value
	a 4-column BED file of the summit position 
Files will use the same path and base name, but change the extension to 
'.bed' and '.summit.bed', respectively.

For gappedPeak files, only one file will be written, a 12-column bed file, 
with extension '.gapped.bed'.

USAGE: peak2bed.pl peak1.narrowPeak [peak1.gappedPeak] ....

OPTIONS
   -i --input <file>    Specify the input file
   -o --output <file>   Specify the output file (default input basename)
   -n --name <text>     Specify the name text (default input basename)
   -s --summit          write a summit bed file (narrowPeak input only)
                           default true, use --nosummit to turn off
   -h --help            print this help
   
   Multiple input and ouput files and names may be specified via options,
   simply repeat as necessary.

DOC

### Options
unless (@ARGV) {
	print $docs;
	exit;
}

my @in_files;
my @out_files;
my @names;
my $write_summit = 1;
my $help;

GetOptions(
	'i|input=s'  => \@in_files,
	'o|output=s' => \@out_files,
	'n|name=s'   => \@names,
	's|summit!'  => \$write_summit,
	'h|help!'    => \$help,
) or die "unrecognized option!\n";

if ($help) {
	print $docs;
	exit;
}

if (@ARGV) {
	push @in_files, @ARGV;
}
unless (@in_files) {
	print "no input files!\n";
	exit 1;
}
if ( @out_files and scalar(@out_files) != scalar(@in_files) ) {
	print "unequal numbers of input files and output files!\n";
	exit 1;
}
if ( @names and scalar(@names) != scalar(@in_files) ) {
	my $n = shift @names;
	@names = ( $n x scalar(@in_files) );
	printf " using '$n' as name for all %d input files\n", scalar(@in_files);
}

# inputs
while ( my $file = shift @in_files ) {

	my ( $basename, $outfile );
	if (@out_files) {
		$outfile = shift @out_files;
	}
	if (@names) {
		$basename = shift @names;
	}

	# Open input narrowPeak
	my $Input = Bio::ToolBox->load_file(
		file  => $file,
		parse => 0,
	);
	unless ($Input) {
		warn "Problem loading file '$file'! skipping\n";
		next;
	}

	# set default values
	unless ($outfile) {
		$outfile = File::Spec->catfile( $Input->path, $Input->basename );
	}
	unless ($basename) {
		$basename = $Input->basename;
	}

	# Sort properly
	$Input->gsort_data if $Input->last_row > 1;

	### NarrowPeak files
	if ( $Input->format eq 'narrowPeak' ) {

		# Open outputs
		if ( $outfile !~ /\.bed(?:\.gz)?$/x ) {
			$outfile .= '.bed';
		}
		my $peak_fh = Bio::ToolBox->write_file($outfile)
			or die "unable to open $outfile for writing! $OS_ERROR\n";
		my ( $summit_file, $summit_fh );
		if ($write_summit) {
			$summit_file = $outfile;
			$summit_file =~ s/bed/summit.bed/;
			$summit_fh = Bio::ToolBox->write_file($summit_file)
				or die "unable to open $summit_file for writing! $OS_ERROR\n";
		}

		# write out
		$Input->iterate(
			sub {
				my $row = shift;

				# write peak
				my $name       = sprintf( "%s.%d", $basename, $row->row_index );
				my $bed_string = $row->bed_string(
					bed  => 5,
					name => $name,
				);
				$peak_fh->printf( "%s\n", $bed_string );

				# write summit if specified
				# easier to do it ourselves here than specify everything
				# to the bed string function
				if ($write_summit) {
					my $start =
						$row->start + $row->value(9);    # start will be 1-based here
					$summit_fh->printf(
						"%s\t%d\t%d\t%s\n", $row->seq_id, $start - 1,
						$start, $name
					);
				}
			}
		);

		# Finish
		$peak_fh->close;
		printf " Wrote %d peaks to $outfile\n", $Input->last_row;
		if ($write_summit) {
			$summit_fh->close;
			print " Wrote summits to $summit_file\n";
		}
	}

	### GappedPeak files
	elsif ( $Input->format eq 'gappedPeak' ) {

		# Open outputs
		if ( $outfile !~ /\.bed(?:\.gz)?$/x ) {
			$outfile .= '.gapped.bed';
		}
		my $peak_fh = Bio::ToolBox->write_file($outfile)
			or die "unable to open $outfile for writing! $OS_ERROR\n";

		# Write out
		$Input->iterate(
			sub {
				my $row  = shift;
				my $name = sprintf( "%s.%d", $basename, $row->row_index );
				my @v    = $row->row_values;
				$peak_fh->printf(
					"%s\n",
					join(
						"\t",  $v[0], $v[1], $v[2], $name, $v[4],
						$v[5], $v[6], $v[7], $v[8], $v[9], $v[10], $v[11]
					)
				);
			}
		);

		# Finish
		$peak_fh->close;
		printf " Wrote %d gapped peaks to $outfile\n", $Input->last_row;
	}

	### BroadPeak files - just in case
	elsif ( $Input->format eq 'broadPeak' ) {

		# Open outputs
		if ( $outfile !~ /\.bed(?:\.gz)?$/x ) {
			$outfile .= '.broad.bed';
		}
		my $peak_fh = Bio::ToolBox->write_file($outfile)
			or die "unable to open $outfile for writing! $OS_ERROR\n";

		# Write out
		$Input->iterate(
			sub {
				my $row  = shift;
				my $name = sprintf( "%s.%d", $basename, $row->row_index );
				my @v    = $row->row_values;
				$peak_fh->printf( "%s\n",
					join( "\t", $v[0], $v[1], $v[2], $name, $v[4] ) );
			}
		);

		# Finish
		$peak_fh->close;
		printf " Wrote %d broad peaks to $outfile\n", $Input->last_row;
	}

	### Empty file
	elsif ( $Input->last_row <= 1 ) {
		warn "File '$file' is empty!\n";
		next;
	}

	### Unrecognized file type
	else {
		warn "File '$file' not a Peak file format!\n";
		next;
	}
}

