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
use Bio::ToolBox::utility qw(format_with_commas);

our $VERSION = 3.3;

# a script to generate mean coverage bedGraph track

my $docs = <<USAGE;
  
A script to generate a chromosomal mean coverage bedGraph track 
to be used in Macs2 as the global control track when there is 
no input for generating a control_lambda chromatin bias track. 

For sparse or low coverage datasets, provide an effective genome 
size to calculate the mean. Otherwise, the non-zero coverage length 
is used, which may be too sparse, resulting in too high of average 
for reliable peak calling of sparse datasets.

It will write out a simple bedGraph representing the genomic mean.

Provide an input file in either bigWig or text bedGraph format. 

VERSION: $VERSION

USAGE: generate_mean_bedGraph.pl -i file.bw -o mean.bdg

OPTIONS:
    -i --in <file>      The input bigWig or bedGraph file
    -o --out <file>     The output bedGraph file (input.global_mean.bdg)
    -g --genome <int>   The effective genome size (bp) to be used
                          (default is actual non-zero coverage)
    -h --help           Display help
        
USAGE

### Options
unless (@ARGV) {
	print $docs;
	exit;
}

my $infile;
my $outfile;
my $genome_size = 0;
my $help;
GetOptions(
	'i|in=s'     => \$infile,
	'o|out=s'    => \$outfile,
	'g|genome=i' => \$genome_size,
	'h|help!'    => \$help,
) or die "unrecognized option! See help\n";

# check options
if ($help) {
	print $docs;
	exit;
}

#### input file
if ( @ARGV and not $infile ) {
	$infile = shift @ARGV;
}
unless ($infile) {
	die "must provide an input bigWig or bedGraph file!\n";
}
my ( $basename, $path, $extension ) = fileparse(
	$infile,
	qw(.bw .bdg .bedgraph .bdg.gz .bedgraph.gz)
);

#### Output file
if ($outfile) {

	# check user provided filename for extension
	unless ( $outfile =~ /\. (?: bdg | bedgraph ) $/xi ) {
		$outfile .= '.bdg';
	}
}
else {
	# generate a new one
	$outfile = $path . $basename . '.global_mean.bdg';
}

#### Process input file	accordingly
if ( $extension eq '.bw' ) {
	process_bw($infile);
}
elsif ( $extension =~ /(?: bdg | bedgraph )/xi ) {
	process_bdg($infile);
}
else {
	die " '$infile' has unrecognized file extension! Try .bw .bdg .bedgraph\n";
}

######## subroutines

sub process_bw {
	my $file = shift;

	# open database
	my $db = Bio::ToolBox->open_database($file)
		or warn " Invalid bigWig file specified! Can't open $file!\n";
	return unless ($db);

	# check type of bigWig
	my $type = ref($db);
	unless ( $type eq 'Bio::DB::Big::File' ) {
		print " Currently not supporting $type adaptors! Only Bio::DB::Big bigWig!\n";
		return;
	}
	unless ( $db->is_big_wig ) {
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
	foreach my $chr ( keys %{$chromhash} ) {

		# iterate over the entire length of the chromosome
		# get 500 intervals per iteration to maybe speed things up?
		my $iter =
			$db->get_intervals_iterator( $chr, 0, $chromhash->{$chr}{length}, 500 );
		while ( my $intervals = $iter->next ) {
			foreach my $interval ( @{$intervals} ) {
				next if $interval->{value} == 0;
				my $length = $interval->{end} - $interval->{start};
				$cov_len   += $length;
				$score_sum += ( $length * $interval->{value} );
			}
		}
	}

	# prepare output
	my $global_mean = calculate_mean( $score_sum, $cov_len );
	my $Data        = Bio::ToolBox->new_data(qw(Chromo Start0 End Score));
	foreach my $c ( sort { $a cmp $b } keys %{$chromhash} ) {

		# the bigWig chromosome hash is random and loses original file order
		# so might as well just sort asciibetically
		$Data->add_row( [ $c, 0, $chromhash->{$c}{length}, $global_mean ] );
	}

	# write out the file
	my $s = $Data->save($outfile);
	print " wrote file $s\n";
}

sub process_bdg {
	my $file = shift;

	# open file handle
	my $infh = Bio::ToolBox->read_file($file)
		or die "unable to open file '$file'!\n";

	# walk through file
	my @chroms;
	my %chromo2length;
	my $current_chrom = q();
	my $cov_len       = 0;
	my $score_sum     = 0;
	while ( my $line = $infh->getline ) {
		chomp $line;
		my @data = split( /\t/, $line );
		next unless scalar(@data) == 4;

		# check chromosome name and sequence
		if ( $data[0] ne $current_chrom ) {
			push @chroms, $data[0];
			$current_chrom = $data[0];
		}

		# always assume that end position is higher than the existing end
		$chromo2length{ $data[0] } = $data[2];

		# ignore zero intervals
		next if $data[3] == 0;

		# add to current data
		my $length = $data[2] - $data[1];
		$cov_len   += $length;
		$score_sum += ( $length * $data[3] );
	}
	$infh->close;
	unless ( $score_sum and $cov_len ) {
		die "no coverage score collected! Was a proper bedGraph provided?\n";
	}

	# prepare output
	my $global_mean = calculate_mean( $score_sum, $cov_len );
	my $Data        = Bio::ToolBox->new_data(qw(Chromo Start0 End Score));
	foreach my $c (@chroms) {
		$Data->add_row( [ $c, 0, $chromo2length{$c}, $global_mean ] );
	}

	# write out the file
	my $s = $Data->save($outfile);
	print " wrote file $s\n";

}

sub calculate_mean {
	my ( $signal, $empirical ) = @_;
	printf " Observed covered genome length is %s bp\n", format_with_commas($empirical);
	my $m = sprintf "%.6f", ( $signal / $empirical );
	print " Empirical genomic mean calculated at $m\n";
	if ($genome_size) {

		# user provided size
		printf " Using provided genome size of %s bp\n", format_with_commas($genome_size);
		$m = sprintf "%.6f", ( $signal / $genome_size );
		print " Using genomic mean calculated at $m\n";
	}
	return $m;
}

