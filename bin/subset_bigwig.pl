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
use File::Which;
use Getopt::Long;
use Parallel::ForkManager;
use File::Path qw(make_path);
use Bio::ToolBox::big_helper qw(generate_chromosome_file);

# variables
my $chrom       = 'chr1';
my $outdir      = './';
my $job         = 2;
my $wig2bw      = sprintf("%s", which 'wigToBigWig');
my $bw2wig      = sprintf("%s", which 'bigWigToWig');
my $help;


### Documentation
my $docs = <<DOC;

A script to subset bigWig files to one specific chromosome. Useful for 
downloading a smaller file to a personal computer for visual evaluation.
New files are written in the specified directory with the same basename 
appended with the chromosome name and extension.

USAGE: subset_bigwig.pl -c chr1 file1.bw file2.bw ...

OPTIONS:
    -c --chrom <text>       The chromosome to subset (default $chrom)
    -o --out <file>         The output directory (default $outdir)
    -j --job <int>          Number of simultaneous jobs, (default $job)
    --wig2bw <path>         ($wig2bw)
    --bw2wig <path>         ($bw2wig)
    --help                  Print documentation

DOC


### Options
unless (@ARGV) {
	print $docs;
	exit;
}
GetOptions(
	'c|chr=s'           => \$chrom,
	'o|out=s'           => \$outdir,
	'wig2bw=s'          => \$wig2bw,
	'bw2wig=s'          => \$bw2wig,
	'j|job=i'           => \$job,
	'help!'             => \$help,
) or die "unrecognized option! See help\n";

# check options
if ($help) {
	print $docs;
	exit;
}


#### bigWig files
my @files = @ARGV;
unless (scalar @files) {
	die "must provide one or more bigWig file!\n";
}


#### Output directory
unless (-e $outdir) {
	make_path($outdir) or die " Cannot make output directory $outdir! $!\n";
	print " Made output directory $outdir\n";
}


#### Chromosome file
my $chromofile = generate_chromosome_file($files[0]) or 
	die "unable to generate a chromsome file from '$files[0]'!\n";



#### Process
my $pm = Parallel::ForkManager->new($job);
foreach my $file (@files) {
	$pm->start and next;
	my (undef, undef, $base) = File::Spec->splitpath($file);
	$base =~ s/\.(?:bw|bigwig)$//i;
	my $outfile = File::Spec->catfile($outdir, $base . ".$chrom.bw");
	my $command = sprintf"%s -chrom=%s %s stdout | %s stdin %s %s", 
		$bw2wig, $chrom, $file, $wig2bw, $chromofile, $outfile;
	if (system($command)) {
		die "something went wrong with command '$command'!!!\n";
	}
	else {
		print " wrote $outfile\n";
	}
	$pm->finish;
}
$pm->wait_all_children;

unlink $chromofile;




