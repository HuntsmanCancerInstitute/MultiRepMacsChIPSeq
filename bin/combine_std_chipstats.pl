#!/usr/bin/perl

use strict;
use IO::File;

unless (@ARGV) {
	print <<END;
A script to combine 
    - Novoalign statistics
    - bam_partial_dedup (or bam_umi_dedup) statistics
    - bam2wig empirical shift determination 
    - Macs2 predicted shift determination

It parses this information from stderr.txt and stdout.txt files from Pysano job
directories. It will write out a single tab-delimited with the numbers. Sample 
names are the directory names.

Usage: $0 <outputfile> 1234X1/ 1234X2/ 1234X3/ ...
END
	exit;
}

my $outfile = shift @ARGV;
my @output;
my %sorter;
my $n = 0;
while (@ARGV) {
	my $dir = shift @ARGV;
	$dir =~ s/\/$//; # strip trailing / if present
	
	# possible stdout values
	my ($total_mapped, $nondup, $dup, $duprate, $keepdup, $trimMeanShift, $extension) = 
		(0, 0, 0, 0, 0, 0, );
	
	# possible stderr values
	my ($novoReads, $novoUnique, $novoUniquePer, $novoMulti, $novoMultiPer, $novoNoMap, 
		$novoNoMapPer, $novoPEmean, $novoPEstdev, $macsFragLength) = 
		(0, 0, '0.0%', 0, '0.0%', 0, '0.0%', 0, 0.0, 0);
	
	# process stdout
	if (-e "$dir/stdout.txt") {
		my $fh = IO::File->new("$dir/stdout.txt");
		while (my $line = $fh->getline) {
			if ($line =~  /^  Total mapped:\s+(\d+)$/) {
				# bam_partial_dedup
				$total_mapped = $1;
			}
			elsif ($line =~ /^\s+(\d+) total mapped alignments$/) {
				# bam_umi_dedup
				$total_mapped = $1;
			}
			elsif ($line =~ /^  Non-duplicate count:\s+(\d+)$/) {
				# bam_partial_dedup
				$nondup = $1;
			}
			elsif ($line =~ /^\s+(\d+) \(\d+\.\d%\) UMI-unique alignments retained$/) {
				# bam_umi_dedup
				$nondup = $1;
			}
			elsif ($line =~ /^  Duplicate count:\s+(\d+)$/) {
				# bam_partial_dedup
				$dup = $1;
			}
			elsif ($line =~ /^\s+(\d+) \(\d+\.\d%\) UMI-duplicate alignments marked$/) {
				# bam_umi_dedup
				$dup = $1;
			}
			elsif ($line =~ /^  Duplication rate:\s+(\d\.\d+)$/) {
				# bam_partial_dedup
				$duprate = $1;
			}
			elsif ($line =~ /^  Retained duplicate count:\s+(\d+)\s*$/) {
				# bam_partial_dedup
				# oops, there may be a space at the end
				$keepdup = $1;
			}
			elsif ($line =~ /^  The trimmed mean shift value is (\d+) /) {
				# bam2wig
				$trimMeanShift = $1;
			}
			elsif ($line =~ /Alignments will be extended by (\d+) bp$/) {
				# bam2wig
				$extension = $1;
			}
		}
		$fh->close;
	}
	
	# process stderr
	if (-e "$dir/stderr.txt") {
		my $fh = IO::File->new("$dir/stderr.txt");
		while (my $line = $fh->getline) {
			if ($line =~  /^#\s+Read Sequences:\s+(\d+)/) {
				$novoReads = $1;
			}
			if ($line =~  /^#\s+Unique Alignment:\s+(\d+) \( ?(\d+\.\d%)\)$/) {
				$novoUnique = $1;
				$novoUniquePer = $2;
			}
			elsif ($line =~ /^#\s+Multi Mapped:\s+(\d+) \( ?(\d+\.\d%)\)$/) {
				$novoMulti = $1;
				$novoMultiPer = $2;
			}
			elsif ($line =~ /^#\s+No Mapping Found:\s+(\d+) \( ?(\d+\.\d%)\)$/) {
				$novoNoMap = $1;
				$novoNoMapPer = $2;
			}
			elsif ($line =~ /^#\s+Mean\s+(\d+),\s+Std Dev\s+(\d+\.\d)$/) {
				$novoPEmean = $1;
				$novoPEstdev = $2;
			}
			elsif ($line =~ /^INFO  @ .+: # predicted fragment length is (\d+) bps $/) {
				$macsFragLength = $1;
			}
		}
		$fh->close;
	}
	
	next unless ($novoReads or $total_mapped);
	if ($nondup and $dup and not $duprate) {
		# this occurs if we're parsing from bam_umi_dedup
		$duprate = sprintf("%.4f", $dup / ($dup + $nondup));
	}
	
	# store
	push @output, join("\t", $dir, $novoReads, $novoUnique, $novoMulti, $novoNoMap, 
		$novoUniquePer, $novoMultiPer, $novoNoMapPer, $novoPEmean, $novoPEstdev, 
		$total_mapped, $nondup, $dup, 
		$duprate, $keepdup, $trimMeanShift, $extension, $macsFragLength);
	
	# sort order
	if ($dir =~ /^(\d+)X(\d+)/) {
		$sorter{num}{$1}{$2} = $n;
	}
	else {
		$sorter{char}{$dir} = $n;
	}
	
	# increment counter
	$n++;
}



# write final
my $fh = IO::File->new($outfile, 'w') or die "can't write to $outfile!\n";
$fh->printf( "%s\n", join("\t", qw(Sample NovoalignTotalReads NovoalignUniqueMapped NovoalignMultiMap 
	NovoalignUnMapped NovoalignUniqueMappedFrac NovoalignMultiMapFrac NovoalignUnMappedFrac 
	NovoalignInsertMean NovoalignInsertStdDev 
	TotalMapped NonDuplicateCount DuplicateCount DuplicateRate RetainedDuplicateCount Bam2WigShift Bam2WigExtension 
	Macs2Extension)));
# sort by experiment ID, requestXsample, e.g. 1234X1
foreach my $i (sort {$a <=> $b} keys %{$sorter{num}}) {
	foreach my $j (sort {$a <=> $b} keys %{$sorter{num}{$i}}) {
		$fh->printf("%s\n", $output[ $sorter{num}{$i}{$j} ]);
	}
}
# sort everything else
foreach my $i (sort {$a cmp $b} keys %{$sorter{char}}) {
	$fh->printf("%s\n", $output[ $sorter{char}{$i} ]);
}
$fh->close;
printf "combined %d samples into $outfile\n", scalar(@output);

