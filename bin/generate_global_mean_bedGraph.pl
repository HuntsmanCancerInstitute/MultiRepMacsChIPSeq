#!/usr/bin/perl

use strict;
use File::Basename qw(fileparse);
use Bio::ToolBox::Data;

# a script to generate mean coverage bedGraph track 

unless (@ARGV) {
print <<USAGE;
  
  A script to generate a globabl mean coverage bedGraph track 
  to be used in Macs2 as the global control track when there is 
  no input for generating a control_lambda chromatin bias track.
  This uses a bigWig file to calculate the global mean and obtain 
  the chromosome sizes. It will write out a simple bedGraph 
  representing the genome with the genomic mean as a single value.
  
  Usage: $0 <file1.bw> ...
  
  It will write out a bedgraph file in the same direcotory and same 
  basename appended with '_global_mean.bdg'.
  
USAGE
exit;
}

foreach my $file (@ARGV) {
	
	# check file
	next unless $file =~ /\.bw$/;
	my ($basename, $path, $extension) = fileparse($file, qq(.bw));
	
	# open bigWig database
	my $bw = Bio::ToolBox::Data->open_database($file);
	
	# prepare output
	my $Data = Bio::ToolBox::Data->new(columns => [qw(Chromo Start0 End Score)]);
	
	# walk through chromosomes
	my $sum = 0;
	my $validCount = 0;
	foreach my $chr ($bw->seq_ids) {
		my $start = 0;
		my $end = $bw->length($chr);
		
		# stats
		my $s = $bw->bigwig->bigWigSummaryArrayExtended($chr,$start,$end,1) or next;
		$s = $s->[0];
		$validCount += $s->{validCount};
		$sum += $s->{sumData};
		
		# add the row to the table, putting in a score of 0 temporarily
		$Data->add_row([$chr, $start, $end, 0]);
	}
	
	# global mean
	my $mean = sprintf("%.4f", $sum/$validCount);
	
	# put in the mean score
	$Data->iterate( sub{ shift->value(3, $mean); } );
	
	# write out the file
	my $global_mean_bdg = $path . $basename . '_expected_mean.bdg';
	my $s = $Data->save($global_mean_bdg);
	print " wrote file $s with global mean coverage of $mean\n";
}

# borrowed from Bio::DB::BigWig
sub _min {
    return $_[0] unless defined $_[1];
    return $_[1] unless defined $_[0];
    return $_[0] < $_[1] ? $_[0] : $_[1];
}
sub _max {
    return $_[0] unless defined $_[1];
    return $_[1] unless defined $_[0];
    return $_[0] < $_[1] ? $_[1] : $_[0];
}




