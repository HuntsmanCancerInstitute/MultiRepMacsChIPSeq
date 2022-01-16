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
use File::Basename qw(fileparse);
use Bio::ToolBox 1.65;

my $VERSION =3;

# a script to generate mean coverage bedGraph track 


my $docs = <<USAGE;
  
  A script to generate a chromosomal mean coverage bedGraph track 
  to be used in Macs2 as the global control track when there is 
  no input for generating a control_lambda chromatin bias track.
  This uses a bigWig or bedGraph coverage file to calculate a global 
  mean. Intervals without coverage are not included in the calculation.
  It will write out a simple bedGraph representing the genome 
  with the respective mean for each chromosome. 

  It will write out a bedgraph file in the same directory with the 
  basename appended with '.global_mean.bdg'. More than one file may 
  be provided at a time.
  
  Version: $VERSION
   
  Usage: $0 <file1.bw> ...
  
USAGE

### Options
unless (@ARGV) {
	print $docs;
	exit;
}

my $infile;
my $outfile;
my $help;
GetOptions(
	'i|in=s'            => \$infile,
	'o|out=s'           => \$outfile,
	'h|help!'           => \$help,
) or die "unrecognized option! See help\n";

# check options
if ($help) {
	print $docs;
	exit;
}


#### input file
if (@ARGV and not $infile) {
	$infile = shift @ARGV;
}
unless ($infile) {
	die "must provide an input bigWig or bedGraph file!\n";
}
my ($basename, $path, $extension) = fileparse($infile, 
	qw(.bw .bdg .bedgraph .bdg.gz .bedgraph.gz));

#### Output file
if ($outfile) {
	# check user provided filename for extension
	unless ($outfile =~ /\.(?:bdg|bedgraph)$/) {
		$outfile .= '.bdg';
	}
}
else {
	# generate a new one
	$outfile = $path . $basename . '.global_mean.bdg';
}



#### Process input file	accordingly
if ($extension eq '.bw') {
	process_bw($infile);
}
elsif ($extension =~ /(?:bdg|bedgraph)/i) {
	process_bdg($infile);
}
else {
	die " '$infile' has unrecognized file extension! Try .bw .bdg .bedgraph\n";
}



######## subroutines

sub process_bw {
	my $file = shift;
	
	# open database
	my $db = Bio::ToolBox->open_database($file) or 
		warn " Invalid bigWig file specified! Can't open $file!\n";
	return unless ($db);
	
	# check type of bigWig
	my $type = ref($db);
	unless ($type eq 'Bio::DB::Big::File') {
		print " Currently not supporting $type adaptors! Only Bio::DB::Big bigWig!\n";
		return;
	}
	unless ($db->is_big_wig) {
		print " database is not a bigWig file!\n";
		return;
	}
	
	# get list of chromosomes
	# this is returned as a hash, so it won't be in same serial order as bigWig file
	# I think that's ok, just maybe not as efficient
	my $chromhash = $db->chroms;
	
	# walk through each chromosome
	my $cov_len   = 0;
	my $score_sum = 0;
	foreach my $chr (keys %$chromhash) {
		# iterate over the entire length of the chromosome
		# get 500 intervals per iteration to maybe speed things up?
		my $iter = $db->get_intervals_iterator($chr, 0, $chromhash->{$chr}{length}, 500);
		while (my $intervals = $iter->next) {
			foreach my $interval (@$intervals) {
				next if $interval->{value} == 0;
				my $length = $interval->{end} - $interval->{start};
				$cov_len += $length;
				$score_sum += ($length * $interval->{value});
			}
		}
	}
	
	# prepare output
	my $Data = Bio::ToolBox->new_data(qw(Chromo Start0 End Score));
	my $global_mean = sprintf "%.6f", ($score_sum / $cov_len);
	foreach my $c (sort {$a cmp $b} keys %$chromhash) {
		# the bigWig chromosome hash is random and loses original file order
		# so might as well just sort asciibetically
		$Data->add_row( [$c, 0, $chromhash->{$c}{length}, $global_mean] );
	}
		
	# write out the file
	my $s = $Data->save($outfile);
	print " wrote file $s\n";
}


sub process_bdg {
	my $file = shift;
	
	# open file handle
	my $infh = Bio::ToolBox->read_file($file) or 
		die "unable to open file '$file'!\n";
	
	# walk through file
	my @chroms;
	my %chromo2length;
	my $current_chrom = '';
	my $cov_len = 0;
	my $score_sum = 0;
	while (my $line = $infh->getline) {
		chomp $line;
		my @data = split("\t", $line);
		next unless scalar(@data) == 4;
		
		# check chromosome name and sequence
		if ($data[0] ne $current_chrom) {
			push @chroms, $data[0];
			$current_chrom = $data[0];
		}
		# always assume that end position is higher than the existing end
		$chromo2length{$data[0]} = $data[2];
		
		# ignore zero intervals
		next if $data[3] == 0;
		
		# add to current data
		my $length = $data[2] - $data[1];
		$cov_len += $length; 
		$score_sum += ($length * $data[3]);
	}
	
	# prepare output
	my $Data = Bio::ToolBox->new_data(qw(Chromo Start0 End Score));
	my $global_mean = sprintf "%.6f", ($score_sum / $cov_len);
	foreach my $c (@chroms) {
		$Data->add_row( [$c, 0, $chromo2length{$c}, $global_mean] );
	}
	
	# write out the file
	my $s = $Data->save($outfile);
	print " wrote file $s\n";
	
}


