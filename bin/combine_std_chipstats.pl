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
use IO::File;
use List::Util qw(sum0);

unless (@ARGV) {
	print <<END;
A script to combine 
    - Novoalign statistics
    - bam_partial_dedup (or bam_umi_dedup) statistics
    - bam2wig empirical shift determination 
    - Macs2 predicted shift determination

It parses this information from standard output and/or error text files from these 
programs. It will also parse the stderr.txt and stdout.txt files from Pysano job
directories. It will write out a single tab-delimited with the numbers. Sample 
names are the given input file names or directory names.

Usage: 
    
    Pysano directories:
        combine_std_chipstats.pl <outputfile> 1234X1/ 1234X2/ ...
    
    Text files:
        combine_std_chipstats.pl <outputfile> file1.txt file2.txt ...
    
END
	exit;
}

# output file
my $outfile = shift @ARGV;
if ($outfile !~ /\.txt$/i) {
	$outfile .= '.txt';
}

# process input files
my @output;
my %sorter;
my $n = 0;
while (@ARGV) {
	
	# look at the given path - it may be a directory or a file
	my $path = shift @ARGV;
	$path =~ s/\/$//; # directory may have a trailing slash - remove if present
	
	# check for files
	my @files;
	if (-d $path) {
		# definitely a directory
		my $stdout = "$path/stdout.txt";
		if (-e $stdout) {
			# a Pysano standard output file
			push @files, $stdout;
		}
		my $stderr = "$path/stderr.txt";
		if (-e $stderr) {
			# a Pysano standard error file
			push @files, $stderr;
		}
	}
	else {
		# assume a readable file
		push @files, $path;
	}
	
	# possible stdout values
	my ($total_mapped, $nondup, $optdup, $optduprate, $workcount, $dup, $duprate, 
		$keepdup, $trimMeanShift, $extension) = 
		(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	
	# possible stderr values
	my ($novoReads, $novoUnique, $novoUniquePer, $novoMulti, $novoMultiPer, $novoNoMap, 
		$novoNoMapPer, $novoPEmean, $novoPEstdev, $macsFragLength) = 
		(0, 0, '0.0%', 0, '0.0%', 0, '0.0%', 0, 0.0, 0);
	
	
	# process files
	foreach my $file (@files) {
		print "  combining $file....\n"; 
		my $fh = IO::File->new($file) or 
			die "unable to read $file! $!\n";
		while (my $line = $fh->getline) {
			
			# total
			if ($line =~  /^  Total mapped:\s+(\d+)$/) {
				# bam_partial_dedup
				$total_mapped = $1;
			}
			elsif ($line =~ /^\s+(\d+) total mapped .*alignments$/) {
				# bam_umi_dedup
				$total_mapped = $1;
			}
			
			# non-duplicate
			elsif ($line =~ /^  Non-duplicate count:\s+(\d+)$/) {
				# bam_partial_dedup
				$nondup = $1;
			}
			elsif ($line =~ /^\s+(\d+) \(\d+\.\d%\) UMI-unique .*alignments retained$/) {
				# bam_umi_dedup
				$nondup = $1;
			}
			
			# optical duplicate
			elsif ($line =~ /^  Optical duplicate count:\s+(\d+)$/) {
				# optical bam_partial_dedup
				$optdup = $1;
			}
			
			# duplicate
			elsif ($line =~ /^  Non-optical duplicate count:\s+(\d+)$/) {
				# non-optical bam_partial_dedup
				$dup = $1;
			}
			elsif ($line =~ /^  Duplicate count:\s+(\d+)$/) {
				# old bam_partial_dedup
				$dup = $1;
			}
			elsif ($line =~ /^\s+(\d+) \(\d+\.\d%\) UMI-duplicate .*alignments/) {
				# bam_umi_dedup
				$dup = $1;
			}
			
			# rate
			elsif ($line =~ /^  Non-optical duplication rate:\s+(\d\.\d+)$/) {
				# non-optical bam_partial_dedup
				$duprate = $1;
			}
			elsif ($line =~ /^  Duplication rate:\s+(\d\.\d+)$/) {
				# regular bam_partial_dedup
				$duprate = $1;
			}
			elsif ($line =~ /^  Optical duplicate rate:\s+(\d\.\d+)$/) {
				# optical bam_partial_dedup
				$optduprate = $1;
			}
			
			# bam partial dedup
			elsif ($line =~ /^  Retained duplicate count:\s+(\d+)\s*$/) {
				# bam_partial_dedup
				# oops, there may be a space at the end
				$keepdup = $1;
			}
			elsif ($line =~ /^  Non-optical working count:\s+(\d+)\s*$/) {
				# non-optical working count bam_partial_dedup
				$workcount = $1;
			}
			
			# bam2wig shift value
			elsif ($line =~ /^  The trimmed mean shift value is (\d+) /) {
				$trimMeanShift = $1;
			}
			elsif ($line =~ /Alignments will be extended by (\d+) bp$/) {
				$extension = $1;
			}
			
			# Novoalign standard error output
			elsif ($line =~  /^#\s+Read Sequences:\s+(\d+)/) {
				$novoReads = $1;
			}
			elsif ($line =~  /^#\s+Unique Alignment:\s+(\d+) \( ?(\d+\.\d)%\)$/) {
				$novoUnique = $1;
				$novoUniquePer = $2;
			}
			elsif ($line =~ /^#\s+Multi Mapped:\s+(\d+) \( ?(\d+\.\d)%\)$/) {
				$novoMulti = $1;
				$novoMultiPer = $2;
			}
			elsif ($line =~ /^#\s+No Mapping Found:\s+(\d+) \( ?(\d+\.\d)%\)$/) {
				$novoNoMap = $1;
				$novoNoMapPer = $2;
			}
			elsif ($line =~ /^#\s+Mean\s+(\d+),\s+Std Dev\s+(\d+\.\d)$/) {
				$novoPEmean = $1;
				$novoPEstdev = $2;
			}
			
			# macs2 predicted fragment length
			elsif ($line =~ /^INFO  @ .+: # predicted fragment length is (\d+) bps $/) {
				$macsFragLength = $1;
			}
		}
		$fh->close;
	}
	
	# infer some counts
	if ($nondup and $dup and not $duprate) {
		# this occurs if we're parsing from bam_umi_dedup
		$duprate = sprintf("%.4f", $dup / ($dup + $nondup));
	}
	if (not $workcount) {
		$workcount = $total_mapped;
	}
	
	# clean up path in case we had a file
	$path =~ s/\.txt$//i;
	
	# store
	push @output, [$path, $novoReads, $novoUnique, $novoMulti, $novoNoMap, 
		$novoUniquePer, $novoMultiPer, $novoNoMapPer, $novoPEmean, $novoPEstdev, 
		$total_mapped, $optdup, $optduprate, $workcount, $nondup, $dup, 
		$duprate, $keepdup, $trimMeanShift, $extension, $macsFragLength];
	
	# sort order
	if ($path =~ /^(\d+)X(\d+)/) {
		# GNomEx identifier
		if (exists $sorter{num}{$1}{$2}) {
			# more than one
			if (ref($sorter{num}{$1}{$2}) eq 'ARRAY') {
				push @{ $sorter{num}{$1}{$2} }, $n;
			}
			else {
				my $first = $sorter{num}{$1}{$2};
				$sorter{num}{$1}{$2} = [$first, $n];
			}
		}
		else {
			$sorter{num}{$1}{$2} = $n;
		}
	}
	else {
		$sorter{char}{$path} = $n;
	}
	
	# increment counter
	$n++;
}

# output headers
my @headers = qw(Sample NovoalignTotalReads NovoalignUniqueMapped NovoalignMultiMap 
	NovoalignUnMapped NovoalignUniqueMappedFrac NovoalignMultiMapFrac NovoalignUnMappedFrac 
	NovoalignInsertMean NovoalignInsertStdDev 
	TotalMapped OpticalDuplicateCount OpticalRate WorkingCount NonDuplicateCount 
	DuplicateCount DuplicateRate RetainedDuplicateCount Bam2WigShift Bam2WigExtension 
	Macs2Extension);

# check for empty columns
	# we skip the sample identifier column, but check everything else
my @columns;
for my $i (1..$#headers) {
	# calculate a sum for this particular column
	# any nonzero sum indicates we have data collected, so this column will be output
	my $check = sum0( map {$output[$_]->[$i]} (0..$#output) );
	push @columns, $i if $check != 0;
}
unshift @columns, 0; # always keep the sample identifier column at beginning

# write final
my $fh = IO::File->new($outfile, 'w') or die "can't write to $outfile!\n";
$fh->printf( "%s\n", join("\t", map { $headers[$_] } @columns) );
# sort by experiment ID, requestXsample, e.g. 1234X1
foreach my $i (sort {$a <=> $b} keys %{$sorter{num}}) {
	foreach my $j (sort {$a <=> $b} keys %{$sorter{num}{$i}}) {
		# check to see if we have an array of multiple samples
		if (ref($sorter{num}{$i}{$j}) eq 'ARRAY') {
			foreach my $n ( @{ $sorter{num}{$i}{$j} }) {
				my $row = $output[$n];
				$fh->printf("%s\n", join("\t", map { $row->[$_] } @columns) )
			}
		}
		else {
			my $row = $output[ $sorter{num}{$i}{$j} ];
			$fh->printf("%s\n", join("\t", map { $row->[$_] } @columns) );
		}
	}
}
# sort everything else
foreach my $i (sort {$a cmp $b} keys %{$sorter{char}}) {
	my $row = $output[ $sorter{char}{$i} ];
	$fh->printf("%s\n", join("\t", map { $row->[$_] } @columns) );
}
$fh->close;
printf "combined %d samples into $outfile\n", scalar(@output);

