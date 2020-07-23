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
use File::Spec;
use Bio::ToolBox 1.65;


### Documentation
my $docs = <<DOC;

A script to convert Macs2 narrowPeak and gappedPeak files into simpler bed files. 
Specifically, discard unused p-value and q-value columns after running Macs2 
bdgpeakcall or bdgbroadcall. 

This script will also resort the peaks properly (numerically, as opposed to 
alphabetically), and rename the peaks using the basename, rather than path.

For narrowPeak files, two files will be written:
	a 5-column BED file of peak intervals, including original score value
	a 4-column BED file of calculated summit position
Files will use the same path and base name, but change the extension to 
'.bed' and '.summit.bed', respectively.

For gappedPeak files, only one file will be written, a 12-column bed file, 
with extension '.gapped.bed'.

No options are available.

USAGE: narrowpeak2bed.pl peak1.narrowPeak [peak1.gappedPeak] ....

DOC




### Options
unless (@ARGV) {
	print $docs;
	exit;
}

# inputs
my @files = @ARGV;
foreach my $f (@files) {
	
	# Open input narrowPeak
	my $Input = Bio::ToolBox->load_file(
		file    => $f,
		parse   => 0,
	);
	unless ($Input) {
		warn "Problem loading file '$f'! skipping\n";
		next;
	}
	my $basename = $Input->basename;
	
	# Sort properly
	$Input->gsort_data;

	### NarrowPeak files
	if ($Input->format eq 'narrowPeak') {
	
		# Open outputs 
		my $peak_file = File::Spec->catfile($Input->path, $basename . ".bed");
		my $peak_fh = Bio::ToolBox->write_file($peak_file) or 
			die "unable to open $peak_file for writing! $!\n";
		my $summit_file = File::Spec->catfile($Input->path, $basename . ".summit.bed");
		my $summit_fh = Bio::ToolBox->write_file($summit_file) or 
			die "unable to open $summit_file for writing! $!\n";
	
		# write out
		$Input->iterate( sub {
			my $row = shift;
		
			# write peak
			my $name = sprintf("%s.%d", $basename, $row->row_index);
			my $bed_string = $row->bed_string(
				bed     => 5,
				name    => $name,
			);
			$peak_fh->printf("%s\n", $bed_string);
		
			# write summit
			# easier to do it ourselves here than specify everything to the bed string function
			my $start = $row->start + $row->value(9); # start will be 1-based here
			$summit_fh->printf("%s\t%d\t%d\t%s\n", $row->seq_id, $start - 1, $start, $name);
		});
	
		# Finish
		$peak_fh->close;
		$summit_fh->close;
		printf " Wrote %d peaks to $peak_file and $summit_file\n", $Input->last_row;
	}
	
	
	### GappedPeak files
	elsif ($Input->format eq 'gappedPeak') {
		
		# Open outputs 
		my $peak_file = File::Spec->catfile($Input->path, $basename . ".gapped.bed");
		my $peak_fh = Bio::ToolBox->write_file($peak_file) or 
			die "unable to open $peak_file for writing! $!\n";
		
		# Write out
		$Input->iterate( sub {
			my $row = shift;
			my $name = sprintf("%s.%d", $basename, $row->row_index);
			my @v = $row->row_values;
			$peak_fh->printf("%s\n", join("\t", $v[0], $v[1], $v[2], $name, $v[4], 
				$v[5], $v[6], $v[7], $v[8], $v[9], $v[10], $v[11]));
		});
		
		# Finish
		$peak_fh->close;
		printf " Wrote %d peaks to $peak_file\n", $Input->last_row;
	}
	
	
	### BroadPeak files - just in case
	elsif ($Input->format eq 'broadPeak') {
		
		# Open outputs 
		my $peak_file = File::Spec->catfile($Input->path, $basename . ".broad.bed");
		my $peak_fh = Bio::ToolBox->write_file($peak_file) or 
			die "unable to open $peak_file for writing! $!\n";
		
		# Write out
		$Input->iterate( sub {
			my $row = shift;
			my $name = sprintf("%s.%d", $basename, $row->row_index);
			my @v = $row->row_values;
			$peak_fh->printf("%s\n", join("\t", $v[0], $v[1], $v[2], $name, $v[4]));
		});
		
		# Finish
		$peak_fh->close;
		printf " Wrote %d peaks to $peak_file\n", $Input->last_row;
	}
	
	
	### Unrecognized file type
	else {
		warn "File '$f' not a Peak file format!\n";
		next;
	}
}



