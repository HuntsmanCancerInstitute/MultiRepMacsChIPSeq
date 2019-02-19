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
use File::Basename qw(fileparse);
use Bio::ToolBox 1.65;
use Bio::ToolBox::db_helper qw(get_chromosome_list);

# a script to generate mean coverage bedGraph track 

unless (@ARGV) {
print <<USAGE;
  
  A script to generate a chromosomal mean coverage bedGraph track 
  to be used in Macs2 as the global control track when there is 
  no input for generating a control_lambda chromatin bias track.
  This uses a bigWig file to calculate each chromosomal mean. It 
  will write out a simple bedGraph representing the genome with 
  the respective mean for each chromosome.
  
  Usage: $0 <file1.bw> ...
  
  It will write out a bedgraph file in the same direcotory and same 
  basename appended with '.global_mean.bdg'.
  
USAGE
exit;
}

foreach my $file (@ARGV) {
	
	# check file
	next unless $file =~ /\.bw$/;
	my ($basename, $path, $extension) = fileparse($file, qq(.bw));
	
	# prepare output
	my $Data = Bio::ToolBox->new_data(qw(Chromo Start0 End));
	
	# walk through chromosomes and fill the table with chromosome name and length
	foreach my $chrsize (get_chromosome_list($file)) {
		# each $chrsize is anonymous array of name and length
		$Data->add_row( [$chrsize->[0], 0, $chrsize->[1]]);
	}
	
	# put in the mean score
	my $i = $Data->add_column('Score');
	$Data->iterate( sub{ 
		my $row = shift;
		my $mean = $row->get_score(dataset => $file, method => 'mean');
		$row->value($i, $mean); 
	} );
	
	# write out the file
	my $out_bdg = $path . $basename . '.global_mean.bdg';
	my $s = $Data->save($out_bdg);
	print " wrote file $s\n";
}




