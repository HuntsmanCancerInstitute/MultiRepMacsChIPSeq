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
use warnings;
use English qw(-no_match_vars);
use Getopt::Long;
use IO::File;
use File::Which;
use File::Spec::Functions qw( curdir splitpath rel2abs );
use Bio::MultiRepChIPSeq::reporter;

our $VERSION = 0.1;

# user variables
my $infile;
my $app = sprintf( "%s", which 'pandoc' );
my $help;


### Documentation
my $docs = <<DOC;

A script to render a markdown report into a self-contained HTML file.
The HTML file will be written in the same directory as the input file.
This uses pandoc, which must be in the environment PATH or explicitly
provided.

VERSION: $VERSION

USAGE: render.pl report.md 

OPTIONS:
  -i --in <file>           The markdown report
  --pandoc <path>          Path to pandoc ($app)
  -h --help                Print documentation
DOC

# Options
unless (@ARGV) {
	print $docs;
	exit;
}
GetOptions(
	'i|in=s'     => \$infile,
	'pandoc=s'   => \$app,
	'h|help!'    => \$help,
) or die "unrecognized option!\n";


# check options
if ($help) {
	print $docs;
	exit;
}
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		print " No input file provided!\n";
		exit 1;
	}
}
unless ($app) {
	print " No pandoc available!\n";
	exit 1;
}

my (undef, $in_path, $md_file) = splitpath($infile);
printf " Preparing to render file '%s' in %s\n", $md_file, $in_path;
my $html_file = $md_file;
$html_file    =~ s/\.md$//;
my $title     = $html_file;
$html_file   .= '.html';



my $current = rel2abs( curdir() );
chdir $in_path;
render();
chdir $current;
exit;


sub render {
	my $html_head  = 'pandoc_head.html';
	my $fh = IO::File->new( $html_head, 'w')
		or die "cannot write extra header file! $OS_ERROR";
	$fh->print( Bio::MultiRepChIPSeq::reporter->pandoc_header );
	$fh->close;
	printf " wrote HTML Head Style file '%s'\n", $html_head;
	
	# pandoc options based on version
	my $cmd_opt2 = sprintf
		"-f gfm -t html --self-contained -H %s -M title=%s -o %s %s",
		$html_head, $title, $html_file, $md_file;
	my $cmd_opt3 = sprintf
		"-f gfm -t html --embed-resources --standalone -H %s -M title=%s -o %s %s",
		$html_head, $title, $html_file, $md_file;
		# these options are for pandoc version 2.x or 3.x
		# which use different options for generating standalone html reports 
		# v3 will respect --self-contained with a warning, but why add extra stress?
	
	# determine version
	my $version = qx($app --version);
	my $command;
	if ( $version =~ /pandoc \s 2 \. (\d+)/x ) {
		if ($1 >= 19) {
			# use version 3 style, deprecates self-contained option
			$command = sprintf "%s %s", $app, $cmd_opt3;
		}
		else {
			$command = sprintf "%s %s", $app, $cmd_opt2;
		}
	}
	elsif ( $version =~ /pandoc \s 3 \. \d+/x ) {
		$command = sprintf "%s %s", $app, $cmd_opt3;
	}
	else {
		# use the older format when unknown????
		$command = sprintf "%s %s", $app, $cmd_opt2;
	}

	# generate html
	print " Executing $command\n";
	system($command);
	if ( -e $html_file ) {
		print " Generated $html_file\n";
		unlink $html_head;
	}
	else {
		print " Something went wrong generating '$html_file'\n";
	}
}
