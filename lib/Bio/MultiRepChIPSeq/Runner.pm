package Bio::MultiRepChIPSeq::Runner;
our $VERSION = 17.8;

=head1 name

Bio::MultiRepChIPSeq::Runner - Object for running the ChIPSeq pipeline

=head1 DESCRIPTION

The main object for running one or more ChIPSeq Experiment Job objects. Jobs 
are added to the Runner, processing commands generated for action on the Job 
input files, and executed by the Runner singly or in parallel. 

=head1 METHODS

=over 4

=item new

This initializes the Runner object, which includes the options hash. 
The options hash should be exported as a reference and provided to 
L<Getopt::Long> for processing of a script command line options. 

=item add_job

Adds a new L<Bio::MultiRepChIPSeq::Job>. Takes the seven parameters 
required therein and passes them directly to respective new() method.

=item list_jobs - returns array of Job objects

=item number_of_jobs - number of Job objects

=item add_command - add a command to list of finished commands

=item progress_file - returns name of progress file

=item check_progress_file - checks and loads contents of progress file if present

=item update_progress_file - updates progress hash and file immediately

=item sample_file - return the filename of the samples replicate/conditions file

=item write_samples_file - write the samples replicate/conditions file

=item run_generate_chr_file - generates the chromosome file

=item execute_commands - executes job commands of external applications

=item _check_command_finished - internal function to check if command finished properly

=item run_input_peak_detection - calls peaks on reference control for exclusion

=item run_dedup - run bam deduplication

=item run_bam_check - checks for deduplicated bam files

=item run_mappable_space_report - calculates mappable space

=item run_bam_fragment_conversion - converts ChIP bam to coverage files

=item run_bam_count_conversion - generates ChIP count bw files

=item run_lambda_control - generates lambda control reference file

=item run_bw_conversion - converts bigWigs to bedGraphs

=item run_bdgcmp - generates fold enrichment and q-value tracks with MACS2

=item run_call_peaks - calls peaks with MACS2

=item run_clean_peaks - convert narrowPeak to bed files

=item run_bdg_conversion - convert bedGraph files to bigWig

=item run_peak_merge - intersect between experiment Job peaks

=item run_rescore - rescore the merged peaks

=item run_efficiency - calculate fragment of reads in peaks

=item run_plot_peaks - run Rscript to plot QC metrics and analysis figures

=item run_cleanup - delete temporary files

=item run_organize - move files into subfolders

=back


=cut

use strict;
use Carp;
use IO::File;
use File::Spec;
use File::Copy;
use File::Path qw(make_path);
use List::Util qw(min);
use Parallel::ForkManager;
use Bio::ToolBox::utility qw(simplify_dataset_name);
use Bio::MultiRepChIPSeq::Job;
use base 'Bio::MultiRepChIPSeq::options';

1;

sub new {
	my $class = shift;
	my $options = $class->init_options();
	
	my $self = {
		opts    			=> $options,
		Jobs    			=> [],
		finished_commands 	=> [],
		pm      			=> undef,
		progress_file 		=> undef,
		sample_file         => undef,
		universal_control   => 0,
	};
	return bless $self, $class;
}

sub version {
	return $VERSION;
}

sub has_universal_control {
	my $self = shift;
	$self->{universal_control} = shift if @_;
	return $self->{universal_control};
}

sub add_job {
	my $self = shift;
	my $Job = Bio::MultiRepChIPSeq::Job->new($self->{opts}, @_);
	if ($Job) {
		push @{ $self->{Jobs} }, $Job;
	}
	else {
		confess "failed to create a Job object!";
	}
	return $Job;
}

sub list_jobs {
	return @{ shift->{Jobs} };
}

sub number_of_jobs {
	return scalar @{ shift->{Jobs} };
}

sub add_command {
	my ($self, $commands) = @_;
	push @{ $self->{finished_commands} }, @$commands;
}

sub progress_file {
	my $self = shift;
	unless (defined $self->{progress_file}) {
		# generate progress file path
		my $pf = File::Spec->catfile($self->dir, $self->out . '.progress.txt');
		$self->{progress_file} = $pf;
		
		# make the progress hash
		$self->{progress} = {
			control_peak  => 0,
			deduplication => 0,
			mappable_size => 0,
			fragment      => 0,
			count         => 0,
			lambda        => 0,
			bw2bdg        => 0,
			bdgcmp        => 0,
			callpeak      => 0,
			cleanpeak     => 0,
			bdg2bw        => 0,
			peakmerge     => 0,
			rescore       => 0,
			efficiency    => 0,
		};
	}
	return $self->{progress_file};
}

sub check_progress_file {
	my $self = shift;
	
	my $pf = $self->progress_file;
	
	if (-e $pf) {
		my $fh = IO::File->new($pf, '<');
		while (my $line = $fh->getline) {
			chomp $line;
			if ($line eq 'bam2wig') {
				# old value!? shouldn't happen but just in case
				$self->{progress}{fragment} = 1;
				$self->{progress}{count} = 1;
			}
			else {
				$self->{progress}{$line} = 1 if exists $self->{progress}{$line};
			}
		}
		$fh->close;
	}
}

sub update_progress_file {
	my $self = shift;
	my $key = shift;
	$self->{progress}{$key} = 1;
	return 1 if $self->dryrun; # just pretend
	
	my $fh = IO::File->new($self->progress_file, '>>') or 
		croak "can't write to progress file! $!\n";
	$fh->print("$key\n");
	$fh->close;
	return 1;
}

sub sample_file {
	my $self = shift;
	if ($self->{sample_file}) {
		return $self->{sample_file};
	}
	else {
		# write the file if we haven't done so yet
		if ($self->number_of_jobs) {
			return $self->write_samples_file;
		}
		else {
			confess "No Jobs to write a sample file!";
		}
	}
}

sub write_samples_file {
	my $self = shift;
	
	# add dataset files
	# we're using count bigWig files for simplicity
	my %name2done;
	my @conditions = ("Replicate\tDataset\n");
	foreach my $Job ($self->list_jobs) {
		foreach my $b ($Job->chip_count_bw) {
			my $name = simplify_dataset_name($b);
			push @conditions, sprintf("%s\t%s\n", $name, $Job->job_name)
		}
		foreach my $b ($Job->control_count_bw) {
			next if exists $name2done{$b};
			my $name = simplify_dataset_name($b);
			push @conditions, "$name\tInput\n";
			$name2done{$b} = 1; # remember it's done
		}
	}
	
	# write samples file
	my $samplefile = File::Spec->catfile($self->dir, $self->out . '_samples.txt');
	unless ($self->dryrun) {
		my $fh = IO::File->new($samplefile, "w");
		foreach (@conditions) {
			$fh->print($_);
		}
		$fh->close;
	}
	
	$self->{sample_file} = $samplefile;
	return $samplefile;
}

sub run_generate_chr_file {
	my $self = shift;
		# this will work regardless if example is bam or bigWig
	print "\n\n======= Generating temporary chromosome file\n";
	my $chromofile = $self->chromofile;
	unless ($chromofile) {
		$chromofile = File::Spec->catfile($self->dir,"chrom_sizes.temp.txt");
		$self->chromofile($chromofile);
	}
	if (-e $chromofile) {
		return $chromofile;
	}
	print " using $chromofile\n";
	unless ($self->printchr_app or $self->dryrun) {
		croak "no print_chromosome_lengths.pl application in path!\n";
	}
	
	# check all source bam and bw files to ensure seqid consistency
	my @sources;
	foreach my $Job ($self->list_jobs) {
		# a job may be just a control job and not have any bam files
		# also possible that no bam files are present
		if (scalar($Job->chip_bams)) {
			push @sources, $Job->chip_bams;
		}
		elsif ($Job->chip_bw and -e $Job->chip_bw) {
				push @sources, $Job->chip_bw;
		}
		if (scalar($Job->control_bams)) {
			push @sources, $Job->control_bams;
		}
	}
	
	# command
	my $command = sprintf("%s --out %s ", 
		$self->printchr_app || 'print_chromosome_lengths.pl', 
		$chromofile
	);
	$command .= sprintf("--chrskip '%s' ", $self->chrskip) if $self->chrskip;
	$command .= join(' ', @sources);
	my $log = $chromofile;
	$log =~ s/\.temp\.txt/.log.txt/;
	$command .= sprintf(" 2>&1 > $log");
	return $self->execute_commands([[$command, $chromofile, $log]]);
}

sub execute_commands {
	my $self = shift;
	my $commands = shift;
	printf "Excecuting %d commands\n", scalar @$commands;
	
		
	# dry run
	if ($self->dryrun) {
		# we just go through the motions here
		foreach my $command (@$commands) {
			printf "=== Job: %s\n", $command->[0];
		}
		$self->add_command($commands);
		return;
	}
	
	# execute jobs
	if ($self->job > 1) {
		
		# get parallel manager
		my $pm;
		if (defined $self->{pm}) {
			$pm = $self->{pm}; # make it reusable
		}
		else {
			# initialize new instance of parallel manager
			$pm = Parallel::ForkManager->new($self->job) or 
				confess "unable to initiate Parallel Forkmanager!";
			$self->{pm} = $pm; # make it reusable
		}
		
		# run commands
		foreach my $command (@$commands) {
			next if $self->_check_command_finished($command, 1);
			printf "=== Job: %s\n", $command->[0];
			
			# check for simple rm commands
			if ($command->[0] =~ /^rm (.+)$/) {
				# we don't need to fork a new process just to execute a rm command
				unlink(split(' ', $1));
				next;
			}
			
			# fork to execute
			$pm->start and next;
			# in child
			system($command->[0]);
			$pm->finish;
		}
		$pm->wait_all_children;
	}
	else {
		foreach my $command (@$commands) {
			next if $self->_check_command_finished($command, 1);
			printf "=== Job: %s\n", $command->[0];
			
			# check for simple rm commands
			if ($command->[0] =~ /^rm (.+)$/) {
				# we don't need to fork a new process just to execute a rm command
				unlink(split(' ', $1));
				next;
			}
			
			# execute
			system($command->[0]);
		}
	}
	
	# check that commands actually produced something
	my @errors;
	foreach my $command (@$commands) {
		unless ( $self->_check_command_finished($command, 0) ) {
			# zero status indicates something went wrong
			push @errors, sprintf "=== ERROR: %s\n", $command->[0];
			if (defined $command->[2]) {
				push @errors, sprintf " See log '%s' for details\n", $command->[2];
			}
		}
	}
	if (@errors) {
		print "\n\n ======= Errors ======\n";
		print " The following jobs did not generate expected output\n";
		foreach (@errors) {
			print $_;
		}
		croak "\nCheck log files for errors\n";
	}
	
	$self->add_command($commands);
}

sub _check_command_finished {
	my $self = shift;
	my ($command, $talk) = @_;
	# returns true if command appears finished
	
	# command bits
	my ($command_string, $command_out, $command_log) = @$command;
	
	# check
	if (length($command_out) and length($command_log)) {
		# both 
		if (-e $command_out and -e $command_log) {
			print "=== Job: $command_string\n    previously finished, have $command_out and $command_log files\n" if $talk;
			return 1;
		}
		elsif (not -e $command_out and -e $command_log) {
			# special instance where there is a log file but no output
			# this can occur with specific applications so check those
			if (index($command_string, $self->bamdedup_app) == 0) {
				# the deduplication command will not write out a bam file if the actual 
				# duplication rate is below the target rate
				# presume this is good
				print "=== Job: $command_string\n    presumed finished, have $command_log file only\n" if $talk;
				return 2;
			}
			elsif (index($command_string, $self->printchr_app) == 0) {
				# the print_chromosome_lengths script will not write output if 
				# sequence orders are not the same
				return 0;
			}
			else {
				# something else? catchall?
				return 0;
			}
		}
		elsif (-e $command_out and not -e $command_log) {
			# we have a output file but not a log file
			print "=== Job: $command_string\n    presumed finished, have $command_out file only, no log\n" if $talk;
			return 3;
		}
	}
	elsif (length($command_out)) {
		if (-e $command_out) {
			print "=== Job: $command_string\n    previously finished, have $command_out\n" if $talk;
			return 4;
		}
	}
	elsif (length($command_out) == 0) {
		# no output files
		if (substr($command_string, 0, 2) eq 'rm') {
			# remove command doesn't leave an output (duh!) or log file
			# gotta check each one
			my $check = 0; 
			foreach my $item (split /\s+/, $command_string) {
				next if $item eq 'rm';
				$check++ if -e $item; # check true if file is present
			}
			if ($check == 0) {
				print "=== Job: $command_string\n    previously finished, target files missing\n" if $talk;
				return 5;
			}
		}
		elsif ($command_string =~ /Rscript/) {
			# plot peaks may leave lots of files or none
			if (-e $command_log) {
				# output log file exists, presume to be finished
				print "=== Job: $command_string\n    presumed finished, have $command_log file\n" if $talk;
				return 2;
			}
		}
	}
	
	# else presume command was not finished
	return 0; 
}

sub run_input_peak_detection {
	my $self = shift;
	return unless ($self->blacklist =~ /^(?:input|control)$/i);
	print "\n\n======= Generating exclusion list from reference control\n";
	if ($self->{progress}{control_peak}) {
		print "\nStep is completed\n";
		return;
	}
	
	# available reference bam files
	my @refbams;
	foreach my $Job ($self->list_jobs) {
		push @refbams, $Job->control_bams;
	}
	my $blacklist;
	if (@refbams) {
		# we have at least one reference bam file to process
		# set the name of the exclusion list file
		$blacklist = File::Spec->catfile($self->dir, sprintf("%s.control_peak", $self->out));
	}
	else {
		# no reference bam files!
		print " No control reference bam files to process. Skipping\n";
		$self->blacklist('');
		$self->update_progress_file('control_peak');
		return;
	}
	
	# check apps
	unless ($self->bam2wig_app =~ /\w+/ or $self->dryrun) {
		croak "no bam2wig.pl application in path!\n";
	}
	unless ($self->macs_app =~ /\w+/ or $self->dryrun) {
		croak "no MACS2 application in path!\n";
	}
	unless ($self->peak2bed_app =~ /\w+/ or $self->dryrun) {
		croak "no peak2bed.pl application in path!\n";
	}
	unless ($self->meanbdg_app =~ /\w+/ or $self->dryrun) {
		croak "no generate_mean_bedGraph.pl application in path!\n";
	}
	
	# generate bam2wig command
	# very little filtering here - we basically want everything
	my $command = sprintf(
		"%s --out %s.bdg --nosecondary --noduplicate --nosupplementary --mean --bdg --bin %s --cpu %s ", 
		$self->bam2wig_app || 'bam2wig.pl', 
		$blacklist,
		$self->chipbin,
		$self->cpu * $self->job, # give it everything we've got, single job
	);
	if ($self->paired) {
		$command .= sprintf("--span --pe --minsize %s --maxsize %s ", 
			$self->minsize, $self->maxsize);
	}
	else {
		$command .= sprintf("--extend --extval %s ", $self->fragsize);
	}
	if ($self->chrskip) {
		$command .= sprintf("--chrskip \'%s\' ", $self->chrskip);
	}
	$command .= sprintf("--in %s ", join(',', @refbams));
	my $logfile = sprintf("%s.out.txt", $blacklist);
	$command .= " 2>&1 > $logfile ";
	
	# add the mean bedgraph file
	$command .= sprintf(" && %s --in %s.bdg --genome %d --out %s.global_mean.bdg 2>&1 >> $logfile ", 
		$self->meanbdg_app || 'generate_mean_bedGraph.pl', 
		$blacklist,
		$self->genome || 0,
		$blacklist
	);
	
	# add the q-value conversion
	$command .= sprintf(" && %s bdgcmp -t %s.bdg -c %s.global_mean.bdg -m qpois -o %s.qvalue.bdg 2>> $logfile ",
		$self->macs_app || 'macs2', 
		$blacklist, 
		$blacklist, 
		$blacklist
	);
	
	# add the peak call
	# we are using hard coded parameters for now, but I think these should be generic enough
	$command .= sprintf(" && %s bdgpeakcall -i %s.qvalue.bdg -c 3 -l 250 -g 500 --no-trackline -o %s.narrowPeak 2>> $logfile",
		$self->macs_app || 'macs2', 
		$blacklist, 
		$blacklist
	);
	
	# execute the call
	$self->execute_commands( [ [$command, "$blacklist.narrowPeak", $logfile], ] );
	
	# check whether peaks were found
	if (-e "$blacklist.narrowPeak" and -s _) {
		# we have identified peaks
		
		$self->blacklist($blacklist . '.bed'); # set the actual output file
		
		# convert to simple bed
		$command = sprintf("%s %s.narrowPeak 2>&1 >> $logfile", 
			$self->peak2bed_app || 'peak2bed.pl', 
			$blacklist
		);
		
		# and clean up
		$command .= sprintf("&& rm %s.bdg %s.global_mean.bdg %s.qvalue.bdg %s.narrowPeak %s.summit.bed",
			$blacklist, 
			$blacklist, 
			$blacklist, 
			$blacklist, 
			$blacklist
		);
		$self->execute_commands( [ [$command, $self->blacklist, $logfile], ] );
	}
	else {
		# no peaks were identified
		print "\n No peaks were found in control\n\n";
		$self->blacklist('');
		
		# clean up
		$command = sprintf("rm %s.bdg %s.global_mean.bdg %s.qvalue.bdg ",
			$blacklist, 
			$blacklist, 
			$blacklist
		);
		$command .= "$blacklist.narrowPeak " if -e "$blacklist.narrowPeak";
		$self->execute_commands( [ [$command, '', ''], ] );
	}
	
	# finished
	$self->update_progress_file('control_peak');
}

sub run_dedup {
	my $self = shift;
	return unless ($self->dedup);
	print "\n\n======= De-duplicating bam files\n";
	if ($self->{progress}{deduplication}) {
		print "\nStep is completed\n";
		return;
	}
	
	### Run de-duplication
	my @commands;
	my %name2done;
	foreach my $Job ($self->list_jobs) {
		push @commands, $Job->generate_dedup_commands(\%name2done);
	}
	if (@commands) {
		$self->execute_commands(\@commands);
		return if $self->dryrun;
	}
	else {
		$self->update_progress_file('deduplication');
		return;
	}
	
	### Collect deduplication statistics
	my @dedupstats;
	push @dedupstats, join("\t", qw(File TotalCount OpticalDuplicateCount DuplicateCount
		NonDuplicateCount DuplicationRate RetainedDuplicateCount));
	foreach my $c (@commands) {
		# initialize counts
		my $total   = 0; # total count
		my $optdup  = 0; # optical duplicate count
		my $nondup  = 0; # non-duplicate count
		my $dup     = 0; # duplicate count
		my $retdup  = 0; # retained duplicate count
		my $duprate = 0; # duplication rate
	
		# open log file and collect stats
		next if (not -e $c->[2]);
		my $fh = IO::File->new($c->[2], 'r');
		while (my $line = $fh->getline) {
			if ($line =~  /^  Total mapped:\s+(\d+)$/) {
				# bam_partial_dedup
				$total = $1;
			}
			elsif ($line =~ /^  Non-duplicate count:\s+(\d+)$/) {
				$nondup = $1;
			}
			elsif ($line =~ /^  Optical duplicate count:\s+(\d+)$/) {
				# optical bam_partial_dedup
				$optdup = $1;
			}
			elsif ($line =~ /^  Non-optical duplicate count:\s+(\d+)$/) {
				# non-optical bam_partial_dedup
				$dup = $1;
			}
			elsif ($line =~ /^  Non-optical duplication rate:\s+(\d\.\d+)$/) {
				# non-optical bam_partial_dedup
				$duprate = $1;
			}
			elsif ($line =~ /^  Retained non-optical duplicate count:\s+(\d+)\s*$/) {
				# bam_partial_dedup
				# oops, there may be a space at the end
				# this might not be found if no deduplication occurred
				$retdup = $1;
			}
		}
		$fh->close;
	
		# name of bam file, extracted from log file
		my (undef, undef, $name) = File::Spec->splitpath($c->[2]);
		$name =~ s/\.dedup\.out\.txt$//;
	
		# store in array
		push @dedupstats, join("\t", $name, $total, $optdup, $dup, $nondup, $duprate, 
			$retdup || $dup);
	}

	# print duplicate stats file
	my $dedupfile = File::Spec->catfile($self->dir, $self->out . '.dedup-stats.txt');
	my $fh = IO::File->new($dedupfile, 'w');
	foreach (@dedupstats) {
		$fh->print("$_\n");
	}
	$fh->close;
	print "\n Wrote deduplication report $dedupfile\n";
	
	$self->update_progress_file('deduplication');
}

sub run_bam_check {
	my $self = shift;
	print "\n\n======= Checking bam files\n";
	foreach my $Job ($self->list_jobs) {
		$Job->find_dedup_bams;
	}
	
	# check for independent peak call and universal control
	if ($self->independent and $self->has_universal_control) {
		my @jobs = $self->list_jobs;
		# universal control is always the first job
		my $control = shift @jobs;
		foreach my $Job (@jobs) {
			$Job->control_use_bams( $control->control_use_bams );
		}
	}
}

sub run_mappable_space_report {
	my $self = shift;
	print "\n\n======= Determining mappable genome space\n";
	
	# first collect all available bam files - can't run without bam files
	my @bamlist;
	foreach my $Job ($self->list_jobs) {
		push @bamlist, ($Job->control_use_bams) || ($Job->control_bams);
		push @bamlist, ($Job->chip_use_bams) || ($Job->chip_bams);
	}
	
	# next get the official genome size from the chromosome file
	my $genome_size = 0;
	if (not $self->dryrun) {
		my $fh = IO::File->new($self->chromofile, 'r');
		while (my $line = $fh->getline) {
			if ($line =~ /\s+(\d+)$/) {
				$genome_size += $1;
			}
		}
		$fh->close;
	}
	
	# check the user supplied value
	if ($self->genome and not $self->dryrun and not $self->{progress}{mappable_size}) {
		
		# double-check the user value
		my $ratio = $self->genome / $genome_size;
		if ($ratio > 1.01) {
			printf "\n User supplied genome size (%d) is larger than actual genome size (%d)!!!\n",
				$self->genome, $genome_size;
			if (@bamlist) {
				printf " Determining actual empirical mappable size\n",
				$self->genome(0);
			}
		}
		elsif ($ratio < 0.6) {
			printf "\n User supplied genome size (%d) is considerably smaller than actual genome size (%d)!\n",
				$self->genome, $genome_size;
			if (@bamlist) {
				printf " Determining actual empirical mappable size\n",
				$self->genome(0);
			}
		}
		else {
			# ratio is somewhere between 50-100% of actual genome size, so assume ok
			printf "\n Using user-specified size of %d\n", $self->genome;
			return;
		}
	}
	elsif ($self->genome and $self->dryrun) {
		printf "\n Pretending that user supplied genome size is ok\n";
		return;
	}
	
	# Check that we have bam files to continue
	unless (@bamlist) {
		# use a default genome size of 90% of actual size
		# this should be close enough for most eukaryotic genomes
		$self->genome(int($genome_size * 0.9));
		printf "\n Input Bam files are not available. Using genome size of %d\n", $self->genome;
		$self->update_progress_file('mappable_size');
		return;
	}
	
	# the output logfile 
	my $logfile = File::Spec->catfile($self->dir, 
		sprintf("%s.mappable.out.txt", $self->out) );
	
	# check if command is finished, otherwise run it
	if ($self->{progress}{mappable_size}) {
		print "\nStep is completed\n";
		# we will read the value below
	}
	else {
		unless ($self->reportmap_app =~ /\w+/ or $self->dryrun) {
			croak "no report_mappable_space.pl application in path!\n";
		}
		my $command = sprintf("%s --cpu %d ", 
			$self->reportmap_app || 'report_mappable_space.pl', 
			$self->cpu * $self->job
		);
			# give the cpu everything we've got, there's only one job
		if ($self->chrskip) {
			$command .= sprintf("--chrskip \'%s\' ", $self->chrskip);
		}
		$command .= join(" ", @bamlist);
		$command .= " 2>&1 > $logfile";
		
		# execute
		# the log file is the output
		$self->execute_commands( [ [$command, $logfile, $logfile] ] );
	}
	
	# Collect results from output file
	# do this regardless whether this was finished previously, since we have to 
	# extract the value into memory anyway
	if (-e $logfile) {
		my $fh = IO::File->new($logfile, 'r') or 
			croak " unable to open mappable report file '$logfile'!";
		while (my $line = $fh->getline) {
			# we're going to use the all mappable space number
			if ($line =~ /All mappable space: ([\d\.]+) Mb/) {
				$self->genome($1 * 1000000);
				last;
			}
		}
		$fh->close;
		if ($self->genome) {
			printf "\n Genome mappable space calculated to be %d bp\n", $self->genome;
		}
		else {
			croak "\n Unable to extract genome mappable space from report '$logfile'!\n";
		}
	}
	elsif ($self->dryrun) {
		print "\n Artificially setting mappable genome size to 100000000 (100 Mb) for dry run purposes\n";
	}
	else {
		croak "\n Genome mappable report log file is not present! Unable to continue!\n";
	}
	
	$self->update_progress_file('mappable_size'); # might be redundant but that's ok
}


sub run_bam_fragment_conversion {
	my $self = shift;
	print "\n\n======= Generating fragment coverage files\n";
	
	
	# Generate commands for each job
	# we need the command regardless of whether it needs to be run or not
	my @commands;
	my %name2done;
	foreach my $Job ($self->list_jobs) {
		push @commands, $Job->generate_bam2wig_frag_commands(\%name2done);
	}
	
	# Execute as necessary
	if ($self->{progress}{fragment}) {
		print "\nStep is completed\n";
		# do NOT return yet
	}
	elsif (@commands) {
		# run programs
		$self->execute_commands(\@commands);
		$self->update_progress_file('fragment');
	}
	
	# skip counting results if dryrun
	if ($self->dryrun) {
		# artificially set target depth
		print "\n Artificially setting target depth to 25 Million for dry run purposes\n";
		$self->targetdep(25);
		return;
	}
	
	# count results
	my %bam2count;
	foreach my $com (@commands) {
		my $log = $com->[2];
		my $fh = IO::File->new($log, 'r') or  # this should open!!!!
			croak "something is wrong! Job completed but unable to open $log!? $!";
		
		my @files; # there may be one or more files processed here
		while (my $line = $fh->getline) {
			chomp $line;
			if ($line =~ /^ Processing files (.+)\.\.\.$/) {
				# the names of files
				@files = split(/, /, $1);
			}
			elsif ($line =~ /^ Normalizing depth based on ([\d,]+) total counted (?:alignments|fragments)$/) {
				# only one file was processed
				# need to grab the name from the list
				my $count = $1;
				unless (exists $bam2count{$files[0]}) {
					$count =~ s/,//g; # remove commas
					$bam2count{$files[0]} = sprintf "%.1f", $count / 1_000_000;
				}
			}
			elsif ($line =~ /^  (.+) had ([\d,]+) total counted (?:alignments|fragments)$/) {
				# multiple files were processed
				my $file = $1;
				my $count = $2;
				unless (exists $bam2count{$file}) {
					$count =~ s/,//g; # remove commas
					$bam2count{$file} = sprintf "%.1f", $count / 1_000_000;
				}
			}
		}
		$fh->close;
	}
	
	# print report
	print "\n Total fragments accepted for analysis\n";
	foreach my $f (sort {$a cmp $b} keys %bam2count) {
		printf "  %6sM  $f\n", $bam2count{$f};
	}
	
	# Calculate minimum target depth to use
	my $targetdep = int( min(values %bam2count) ) || 1; # just in case!
	if (defined $self->targetdep) {
		printf "\n WARNING!!! Calculated target sequence depth of %d Million is overridden by manually set value %d\n",
			$targetdep, $self->targetdep;
	}
	else {
		$self->targetdep($targetdep);
		printf "\n Setting target sequence depth to %d Million\n", $self->targetdep;
	}
}

sub run_bam_count_conversion {
	my $self = shift;
	print "\n\n======= Generating fragment count files\n";
	if ($self->{progress}{count}) {
		print "\nStep is completed\n";
		return;
	}
	
	# generate commands and run
	my @commands;
	my %name2done;
	foreach my $Job ($self->list_jobs) {
		push @commands, $Job->generate_bam2wig_count_commands(\%name2done);
	}
	if (@commands) {
		# run programs
		$self->execute_commands(\@commands);
	}
	
	$self->update_progress_file('count');
}

sub run_lambda_control {
	my $self = shift;
	print "\n\n======= Generating control files\n";
	if ($self->{progress}{lambda}) {
		print "\nStep is completed\n";
		return;
	}
	
	# generate commands
	my @commands;
	my %name2done;
	my @jobs = $self->list_jobs;
	if ($self->has_universal_control) {
		# only need to process the first job - the universal control
		my $con_Job = shift @jobs;
		push @commands, $con_Job->generate_lambda_control_commands(\%name2done);
	}
	else {
		foreach my $Job (@jobs) {
			push @commands, $Job->generate_lambda_control_commands(\%name2done);
		}
	}
	if (@commands) {
		$self->execute_commands(\@commands);
	}
	$self->update_progress_file('lambda');
}

sub run_bw_conversion {
	my $self = shift;
	print "\n\n======= Converting Fragment bigWig files to bedGraph\n";
	if ($self->{progress}{bw2bdg}) {
		print "\nStep is completed\n";
		return;
	}
	my @commands;
	my %name2done;
	foreach my $Job ($self->list_jobs) {
		push @commands, $Job->convert_bw_to_bdg(\%name2done);
	}
	if (@commands) {
		$self->execute_commands(\@commands);
	}
	$self->update_progress_file('bw2bdg');
}

sub run_bdgcmp {
	my $self = shift;
	print "\n\n======= Generate enrichment files\n";
	if ($self->{progress}{bdgcmp}) {
		print "\nStep is completed\n";
		return;
	}
	
	# list of jobs to work on
	my @jobs = $self->list_jobs;
	if ($self->has_universal_control) {
		# remove the control job
		shift @jobs;
	}
	
	# generate commands
	my @commands;
	foreach my $Job (@jobs) {
		push @commands, $Job->generate_enrichment_commands();
	}
	$self->execute_commands(\@commands) if @commands;
	$self->update_progress_file('bdgcmp');
}

sub run_call_peaks {
	my $self = shift;
	print "\n\n======= Call peaks\n";
	if ($self->{progress}{callpeak}) {
		print "\nStep is completed\n";
		return;
	}
	
	# list of jobs to work on
	my @jobs = $self->list_jobs;
	if ($self->has_universal_control) {
		# remove the control job
		shift @jobs;
	}
	
	# generate commands based on whether we're doing joint or independent calls
	my @commands;
	foreach my $Job (@jobs) {
		if ($self->independent) {
			push @commands, $Job->generate_independent_peakcall_commands;
		}
		else {
			push @commands, $Job->generate_peakcall_commands;
		}
	}
	$self->execute_commands(\@commands) if @commands;
	
	# independent calls must be merged
	if ($self->independent) {
		@commands = ();
		foreach my $Job (@jobs) {
			push @commands, $Job->generate_independent_merge_peak_commands;
		}
		if (@commands) {
			$self->execute_commands(\@commands);
		}
	}
	
	$self->update_progress_file('callpeak');
}

sub run_clean_peaks {
	my $self = shift;
	return if $self->independent; # not necessary here
	print "\n\n======= Cleaning peak files\n";
	if ($self->{progress}{cleanpeak}) {
		print "\nStep is completed\n";
		return;
	}
	
	# list of jobs to work on
	my @jobs = $self->list_jobs;
	if ($self->has_universal_control) {
		# remove the control job
		shift @jobs;
	}
	
	my @commands;
	foreach my $Job (@jobs) {
		push @commands, $Job->generate_cleanpeak_commands;
	}
	$self->execute_commands(\@commands);
	$self->update_progress_file('cleanpeak');
}

sub run_bdg_conversion {
	my $self = shift;
	print "\n\n======= Converting bedGraph files\n";
	if ($self->{progress}{bdg2bw}) {
		print "\nStep is completed\n";
		return;
	}
	my @commands;
	my %name2done;
	foreach my $Job ($self->list_jobs) {
		push @commands, $Job->generate_bdg2bw_commands(\%name2done);
	}
	$self->execute_commands(\@commands);
	$self->update_progress_file('bdg2bw');
}

sub run_peak_merge {
	my $self = shift;
	return if $self->number_of_jobs == 1; # no sense merging one job!
	print "\n\n======= Merging called Peak files\n";
	if ($self->{progress}{peakmerge}) {
		print "\nStep is completed\n";
		return;
	}
	unless ($self->bedtools_app =~ /\w+/ or $self->dryrun) {
		croak "no bedtools application in path!\n";
	}
	unless ($self->intersect_app =~ /\w+/ or $self->dryrun) {
		croak "no intersect_peaks.pl application in path!\n";
	}
	
	# list of jobs to work on
	my @jobs = $self->list_jobs;
	if ($self->has_universal_control) {
		# remove the control job
		shift @jobs;
	}
	
	my @commands;
	
	# narrowPeaks
	my $merge_file = File::Spec->catfile($self->dir, $self->out);
	my $command = sprintf("%s --name %s_merge --bed %s --out %s --genome %s ", 
		$self->intersect_app || 'intersect_peaks.pl', 
		$self->out, 
		$self->bedtools_app || 'bedtools', 
		$merge_file, 
		$self->chromofile
	);
	my $command_check = length($command);
	my $count_check = 0;
	foreach my $Job (@jobs) {
		if (
			$Job->clean_peak and 
			( (-e $Job->clean_peak and -s _ > 0) or $self->dryrun )
		) {
			# only add the file if it exists and non-zero in length
			# or just fake it if we're running a dry run
			$command .= sprintf("%s ", $Job->clean_peak);
			$count_check++;
		}
	}
	if (length($command) > $command_check and $count_check > 1) {
		my $log = $merge_file . '.merge.out.txt';
		$command .= sprintf("2>&1 > %s ", $log);
		push @commands, [$command, $merge_file . '.bed', $log]; 
			# this will have multiple outputs, but one is just a .bed file
	}
	else {
		unless ($self->dryrun) {
			print "One or fewer narrow peak files found, nothing to merge!\n";
		}
	}
	
	# broadPeaks
	if ($self->broad) {
		my $merge2_file = File::Spec->catfile($self->dir, $self->out . "_broad");
		my $command2 = sprintf("%s --name %s_gapmerge --bed %s --out %s --genome %s ", 
			$self->intersect_app || 'intersect_peaks.pl', 
			$self->out, 
			$self->bedtools_app || 'bedtools', 
			$merge2_file, 
			$self->chromofile
		);
		my $command2_check = length($command2);
		my $count2_check = 0;
		
		foreach my $Job (@jobs) {
			if (
				$Job->clean_gappeak and 
				( (-e $Job->clean_gappeak and -s _  > 0) or $self->dryrun )
			) {
				# only add the file if it exists and non-zero in length
				# or just fake it if we're running a dry run
				$command2 .= sprintf("%s ", $Job->clean_gappeak);
				$count2_check++;
			}
		}
		if (length($command2) > $command2_check and $count2_check > 1) {
			my $log2 = $merge2_file . '.merge.out.txt';
			$command2 .= sprintf("2>&1 > %s ", $log2);
			push @commands, [$command2, $merge2_file . '.bed', $log2];
				# this will have multiple outputs, but one is just a .bed file
		}
		else {
			unless ($self->dryrun) {
				print "One or fewer gapped peak files found, nothing to merge!\n";
			}
		}
	}
	
	if (@commands) {
		$self->execute_commands(\@commands);
		$self->update_progress_file('peakmerge');
	}
}

sub run_rescore {
	my $self = shift;
	print "\n\n======= Re-scoring all peaks\n";
	if ($self->{progress}{rescore}) {
		print "\nStep is completed\n";
		return;
	}
	unless ($self->getdata_app =~ /\w+/ or $self->dryrun) {
		croak "no get_datasets.pl application in path!\n";
	}
	unless ($self->getrel_app =~ /\w+/ or $self->dryrun) {
		croak "no get_relative_data.pl application in path!\n";
	}
	
	# prepare filenames
	my $input;
	if ($self->number_of_jobs > 1) {
		# multiple jobs, regenerate the merged peak file name
		$input = File::Spec->catfile($self->dir, $self->out . '.bed');
		if (not $self->dryrun and not -e $input) {
			print "No merged peak file '$input'!\n";
			return;
		}
	}
	else {
		# only one job, so take its peak file
		$input = ($self->list_jobs)[0]->clean_peak;
		if (not $self->dryrun and not -e $input) {
			print "No peak file '$input'!\n";
			return;
		}
	}
	my $output1 = File::Spec->catfile($self->dir, $self->out . '_meanQvalue.txt');
	my $output2 = File::Spec->catfile($self->dir, $self->out . '_meanLog2FE.txt');
	my $output3 = File::Spec->catfile($self->dir, $self->out . '_counts.txt');
	my $output4 = File::Spec->catfile($self->dir, $self->out . '_profile_fragment.txt');
	my $output5 = File::Spec->catfile($self->dir, $self->out . '_profile_log2FE.txt');
	my $output6 = File::Spec->catfile($self->dir, $self->out . '_genome_counts.txt.gz');
	
	# generate four get_dataset and two get_relative commands
	# go ahead and make the fourth genome-wide command, even though we may not use it
	my @command_lengths;
	my $command1 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
		$self->getdata_app || 'get_datasets.pl', 
		$self->cpu, 
		$input, 
		$output1
	);
	push @command_lengths, length($command1);
	my $command2 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
		$self->getdata_app || 'get_datasets.pl', 
		$self->cpu, 
		$input, 
		$output2
	);
	push @command_lengths, length($command2);
	my $command3 = sprintf("%s --method sum --cpu %s --in %s --out %s --format 0 ",
		$self->getdata_app || 'get_datasets.pl', 
		$self->cpu, 
		$input, 
		$output3
	);
	push @command_lengths, length($command3);
	my $command4 = sprintf("%s --method mean --cpu %s --in %s --out %s --win %s --num 25 --pos m --long --format 3 --groups --sum ",
		$self->getrel_app || 'get_relative_data.pl', 
		$self->cpu, 
		$input, 
		$output4, 
		$self->binsize
	);
	push @command_lengths, length($command4);
	my $command5 = sprintf("%s --method mean --cpu %s --in %s --out %s --win %s --num 25 --pos m --long --format 3 --groups --sum ",
		$self->getrel_app || 'get_relative_data.pl', 
		$self->cpu, 
		$input, 
		$output5, 
		$self->binsize
	);
	push @command_lengths, length($command5);
	my $command6 = sprintf("%s --method sum --cpu %s --feature genome --win %d --discard %s --out %s --format 0 ",
		$self->getdata_app || 'get_datasets.pl', 
		$self->cpu, 
		$self->genomewin, 
		$self->discard, 
		$output6
	);
	push @command_lengths, length($command6);
	
	# add dataset files
	my %name2done;
	foreach my $Job ($self->list_jobs) {
		if ($Job->qvalue_bw) {
			$command1 .= sprintf("--data %s ", $Job->qvalue_bw);
		}
		if ($Job->logfe_bw) {
			$command2 .= sprintf("--data %s ", $Job->logfe_bw);
			$command5 .= sprintf("--data %s ", $Job->logfe_bw);
		}
		if ($Job->chip_bw) {
			$command4 .= sprintf("--data %s ", $Job->chip_bw);
		}
		foreach my $b ($Job->chip_count_bw) {
			next if exists $name2done{$b};
			$command3 .=  "--data $b ";
			$command6 .=  "--data $b ";
			$name2done{$b} = 1;
		}
		foreach my $b ($Job->control_count_bw) {
			next if exists $name2done{$b};
			$command3 .=  "--data $b ";
			$command6 .=  "--data $b ";
			$name2done{$b} = 1;
		}
	}
	
	
	# add log outputs to commands
	# check the length of each command string to make sure it's valid
	my @commands;
	
	if (length $command1 > shift @command_lengths) {
		my $log = $output1;
		$log =~ s/txt$/out.txt/;
		$command1 .= " 2>&1 > $log";
		push @commands, [$command1, $output1, $log];
	}
	if (length $command2 > shift @command_lengths) {
		my $log = $output2;
		$log =~ s/txt$/out.txt/;
		$command2 .= " 2>&1 > $log";
		push @commands, [$command2, $output2, $log];
	}
	if (length $command3 > shift @command_lengths) {
		my $log = $output3;
		$log =~ s/txt$/out.txt/;
		$command3 .= " 2>&1 > $log";
		push @commands, [$command3, $output3, $log];
	}
	if (length $command4 > shift @command_lengths) {
		my $log = $output4;
		$log =~ s/txt$/out.txt/;
		$command4 .= " 2>&1 > $log";
		push @commands, [$command4, $output4, $log];
	}
	if (length $command5 > shift @command_lengths) {
		my $log = $output5;
		$log =~ s/txt$/out.txt/;
		$command5 .= " 2>&1 > $log";
		push @commands, [$command5, $output5, $log];
	}
	if (length $command6 > shift @command_lengths) {
		if ($self->genomewin) {
			# user has given an actual genome window size, so we'll run this command
			my $log = $output6;
			$log =~ s/txt\.gz$/out.txt/;
			$command6 .= " 2>&1 > $log";
			push @commands, [$command6, $output6, $log];
		}
	}
	
	
	# broad peak rescore
	my ($output7, $output8, $output9);
	if ($self->broad) {
		# generate broad file names
		my $input2;
		if (scalar(@{$self->{Jobs}}) > 1) {
			# merged broad file name
			$input2 = File::Spec->catfile($self->dir, $self->out . '_broad.bed');
			if (not $self->dryrun and not -e $input2) {
				croak "unable to find merged gapped Peak bed file '$input2'!\n";
			}
		}
		else {
			# just one job, take its broad peak file
			$input2 = $self->{Jobs}->[0]->{clean_gappeak};
			if (not $self->dryrun and not -e $input2) {
				croak "unable to find gapped Peak bed file '$input2'!\n";
			}
		}
		$output7 = File::Spec->catfile($self->dir, $self->out . '_broad_meanQvalue.txt');
		$output8 = File::Spec->catfile($self->dir, $self->out . '_broad_meanLog2FE.txt');
		$output9 = File::Spec->catfile($self->dir, $self->out . '_broad_counts.txt');
		
		# generate three get_dataset commands
		# we won't run get_relative data for gapped peaks
		my $command7 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
			$self->getdata_app || 'get_datasets.pl', 
			$self->cpu, 
			$input2, 
			$output7
		);
		push @command_lengths, length($command7);
		my $command8 = sprintf("%s --method mean --cpu %s --in %s --out %s --format 3 ",
			$self->getdata_app || 'get_datasets.pl', 
			$self->cpu, 
			$input2, 
			$output8
		);
		push @command_lengths, length($command8);
		my $command9 = sprintf("%s --method sum --cpu %s --in %s --out %s --format 0 ",
			$self->getdata_app || 'get_datasets.pl', 
			$self->cpu, 
			$input2, 
			$output9
		);
		push @command_lengths, length($command8);
		
		# add dataset files for broad peaks
		%name2done = ();
		foreach my $Job ($self->list_jobs) {
			if ($Job->qvalue_bw) {
				$command7 .= sprintf("--data %s ", $Job->qvalue_bw);
			}
			if ($Job->logfe_bw) {
				$command8 .= sprintf("--data %s ", $Job->logfe_bw);
			}
			foreach my $b ($Job->chip_count_bw) {
				$command9 .=  "--data $b ";
			}
			foreach my $b ($Job->control_count_bw) {
				next if exists $name2done{$b};
				$command9 .=  "--data $b ";
				$name2done{$b} = 1; # remember it's done
			}
		}
	
		# add log outputs to commands
		if (length($command7) > shift @command_lengths) {
			my $log = $output7;
			$log =~ s/txt$/out.txt/;
			$command7 .= " 2>&1 > $log";
			push @commands, [$command7, $output7, $log];
		}
		if (length($command8) > shift @command_lengths) {
			my $log = $output8;
			$log =~ s/txt$/out.txt/;
			$command8 .= " 2>&1 > $log";
			push @commands, [$command8, $output8, $log];
		}
		if (length($command9) > shift @command_lengths) {
			my $log = $output9;
			$log =~ s/txt$/out.txt/;
			$command9 .= " 2>&1 > $log";
			push @commands, [$command9, $output9, $log];
		}
	}
	
	
	
	### Execute data collection commands
	$self->execute_commands(\@commands);
	
	
	
	
	### Replicate Merge
	# do this here so that we still know the output commands
	# must be done after the data collection anyway, so can't be merged above
	if ($self->repmean) {
		# we will generate count means of the replicates
		unless ($self->combrep_app =~ /\w+/ or $self->dryrun) {
			croak "no combine_replicate_data.pl application in path!\n";
		}
		print "\n\n======= Generating sample replicate means\n";
		my @commands2;
		
		# narrowPeak counts
		my $output3m = File::Spec->catfile($self->dir, $self->out . '_meanCounts.txt');
		my $command10 = sprintf("%s --in %s --out %s --sample %s --method mean --format 0 ", 
			$self->combrep_app || 'combine_replicate_data.pl', 
			$output3, 
			$output3m, 
			$self->sample_file
		);
		my $log = $output3m;
		$log =~ s/txt/out.txt/;
		$command10 .= " 2>&1 > $log";
		push @commands2, [$command10, $output3m, $log];
		
		# genome counts
		if ($self->genomewin) {
			my $output6m = File::Spec->catfile($self->dir, $self->out . 
				'_genome_meanCounts.txt.gz');
			my $command11 = sprintf("%s --in %s --out %s --sample %s --method mean --format 0 ", 
				$self->combrep_app || 'combine_replicate_data.pl', 
				$output6, 
				$output6m, 
				$self->sample_file
			);
			$log = $output6m;
			$log =~ s/txt\.gz/out.txt/;
			$command11 .= " 2>&1 > $log";
			push @commands2, [$command11, $output6m, $log];
		}
		
		# broadPeak Counts
		if ($self->broad) {
			my $output9m = File::Spec->catfile($self->dir, $self->out . 
				'_broad_meanCounts.txt');
			my $command12 = sprintf("%s --in %s --out %s --sample %s --method mean --format 0 ", 
				$self->combrep_app || 'combine_replicate_data.pl', 
				$output9, 
				$output9m, 
				$self->sample_file
			);
			$log = $output9m;
			$log =~ s/txt/out.txt/;
			$command12 .= " 2>&1 > $log";
			push @commands2, [$command12, $output9m, $log];
		}
		
		# execute
		$self->execute_commands(\@commands2);
	}
	
	$self->update_progress_file('rescore');
}

sub run_efficiency {
	my $self = shift;
	print "\n\n======= Scoring peak calls for ChIP efficiency\n";
	if ($self->{progress}{efficiency}) {
		print "\nStep is completed\n";
		return;
	}
	unless ($self->geteff_app =~ /\w+/ or $self->dryrun) {
		croak "no get_chip_efficiency.pl application in path!\n";
	}
	
	### ChIP efficiency
	my @commands;
	
	# universal control counts
	# I have to search for these explicitly, since they're not associated with ChIP job
	my @universal_counts;
	foreach my $Job ($self->list_jobs) {
		next unless (
			not defined $Job->clean_peak and 
			scalar($Job->control_count_bw) > 0
		);
		push @universal_counts, ($Job->control_count_bw);
	}
	
	# collect count files for each job
	foreach my $Job ($self->list_jobs) {
		next unless defined $Job->clean_peak;
		next unless scalar($Job->chip_count_bw) > 0;
		my $output = File::Spec->catfile($self->dir, $Job->job_name . '.efficiency.txt');
		my $command = sprintf("%s --in %s --group %s --out %s --cpu %d ", 
			$self->geteff_app || 'get_chip_efficiency.pl', 
			$Job->clean_peak, 
			$self->sample_file, 
			$output, 
			$self->cpu
		);
		# add count files, we should have at least one chip and one control
		foreach my $b ($Job->chip_count_bw) {
			$command .= "$b ";
		}
		foreach my $b ($Job->control_count_bw) {
			$command .= "$b ";
		}
		foreach my $b (@universal_counts) {
			$command .= "$b ";
		}
		my $log = $output;
		$log =~ s/txt/out.txt/;
		$command .= sprintf("2>&1 > $log");
		push @commands, [$command, $output, $log];
	}
	
	# execute the efficiency commands if present
	unless (@commands) {
		# nothing to do? so return
		$self->update_progress_file('efficiency');
		return;
	}
	$self->execute_commands(\@commands);
	
	# proceed no further if dry run
	return if $self->dryrun;
	
	# merge the efficiency outputs into one
	if (scalar @commands > 1) {
		# we're doing this manually to capture all the metadata lines
		my @combined_eff_data;
		my @combined_eff_meta;
		my $eff_header;
		foreach my $c (@commands) {
			my $f = $c->[1];
			unless (-e $f) {
				print " Missing efficiency file $f! Skipping\n";
				next;
			}
			my $fh = IO::File->new($f, 'r') or next;
			while (my $line = $fh->getline) {
				if (substr($line, 0, 1) eq '#') {
					push @combined_eff_meta, $line;
				}
				elsif ($line =~ /^Replicate/) {
					$eff_header = $line unless $eff_header;
				}
				else {
					push @combined_eff_data, $line;
				}
			}
			$fh->close;
			unlink $c->[1];
		}
		
		# write merged file
		if (@combined_eff_data) {
			my $combined_eff_out = File::Spec->catfile($self->dir, $self->out . '.chip_efficiency.txt');
			my $fh = IO::File->new($combined_eff_out, 'w');
			foreach (@combined_eff_meta) {
				$fh->print($_);
			}
			$fh->print($eff_header);
			foreach (@combined_eff_data) {
				$fh->print($_);
			}
			$fh->close;
			print "\nWrote combined ChIP efficiency file $combined_eff_out\n";
		}
		else {
			print " No efficiency files combined!\n";
			return;
		}
	}
	elsif (scalar @commands == 1) {
		# just one, so rename it for consistency sake
		my $eff_out = File::Spec->catfile($self->dir, $self->out . '.chip_efficiency.txt');
		move($commands[0]->[1], $eff_out);
	}
	
	$self->update_progress_file('efficiency');
}

sub run_plot_peaks {
	my $self = shift;
	return unless ($self->plot);
	print "\n\n======= Plotting Peak figures\n";
	if ($self->rscript_app !~ /\w+/ or $self->plotpeak_app !~ /\w+/) {
		# don't die here, just return safely - R is hard
		print "Rscript or plot_peak_figures.R script not defined!\n";
		return;
	}
	my $outbase = File::Spec->catfile($self->dir, $self->out);
	my $command = sprintf("%s --verbose %s --input %s ", $self->rscript_app, 
		$self->plotpeak_app, $outbase);
	my $log = $outbase . '.plot_figures.out.txt';
	$command .= " > $log 2>&1"; # error redirect needs to last with Rscript
	
	# there are multiple possible output files, including none, from this script
	# hard to predict, so don't put any
	my @commands = ( [$command, '', $log] );
	
	# broad peak
	if ($self->broad) {
		my $outbase2 = File::Spec->catfile($self->dir, $self->out) . '_broad';
		my $command2 = sprintf("%s --verbose %s --input %s ", $self->rscript_app, 
			$self->plotpeak_app, $outbase2);
		my $log2 = $outbase2 . '.plot_figures.out.txt';
		$command2 .= " > $log2 2>&1";
		push @commands, [$command2, '', $log2];
	}
	
	# add independent peak calls
	if ($self->independent) {
		foreach my $Job ($self->list_jobs) {
			next if (scalar($Job->rep_peaks) == 1); # nothing to compare
			my $jobbase = File::Spec->catfile($self->dir, $Job->job_name);
			my $command = sprintf("%s --verbose %s --input %s ", $self->rscript_app, 
				$self->plotpeak_app, $jobbase);
			my $example = $jobbase . '.jaccard.png';
			my $log = $jobbase . '.plot_figures.out.txt';
			$command .= " > $log 2>&1"; # error redirect needs to last with Rscript
			push @commands, [$command, $example, $log];
			
			# broad peak
			if ($self->broad) {
				my $jobbase2 = $Job->clean_gappeak;
				$jobbase2 =~ s/\.bed$//;
				my $command2 = sprintf("%s --verbose %s --input %s ", $self->rscript_app, 
					$self->plotpeak_app, $jobbase2);
				$example = $jobbase2 . '.jaccard.png';
				my $log2 = $jobbase2 . '.plot_figures.out.txt';
				$command2 .= " > $log2 2>&1";
				push @commands, [$command2, $example, $log2];
			}
	
		}
	}
	
	$self->execute_commands(\@commands);
}

sub print_config {
	my $self = shift;
	my $capture = shift || 0;
	my @output;
	
	# Samples
	push @output, "\n\n======= Samples\n";
	my @names = $self->name;
	for my $i (0 .. $#names) {
		push @output, sprintf " %s ChIP: %s\n", $names[$i], ($self->chip)[$i]; 
		push @output, sprintf " %s Control: %s\n", $names[$i], ($self->control)[$i] || q();
	}
	
	# Run parameters
	push @output, "\n\n======= Configuration\n";
	foreach my $k (sort {$a cmp $b} keys %{$self->{opts}} ) {
		if (ref($self->{opts}->{$k}) eq 'ARRAY') {
			# next if $k =~ /^(?:chip|control|name)$/;
			push @output, sprintf "%12s  %s\n", $k, 
				join(",\n              ", @{$self->{opts}->{$k}} );
		}
		else {
			push @output, sprintf "%12s  %s\n", $k, $self->{opts}->{$k};
		}
	}
	
	# Finish
	if ($capture) {
		return @output;
	}
	else {
		print @output;
	}
	return;
}


sub run_cleanup {
	my $self = shift;
	return if $self->dryrun;
	print "\n\n======= Combining log files\n";
	my @output;
	
	# Version
	push @output, "======== ChIPSeq multi-replicate pipeline ==========\n";
	push @output, "\nProgram $0\n";
	push @output, "Version $VERSION\n";
	
	# Configuration
	push @output, $self->print_config(1); # capture it
	
	# Combine known output logs
	push @output, "\n\n======= Job Logs\n";
	foreach my $command (@{ $self->{finished_commands} }) {
		# job command string
		push @output, sprintf "\n\n=== Job: %s\n", $command->[0];
		
		# job log file
		if ($command->[2] and -e $command->[2]) {
			if (-s _ ) {
				# file is not empty
				my $fh = IO::File->new($command->[2], 'r') or next;
				push @output, <$fh>;
				$fh->close;
			}
			unlink $command->[2];
		}
	}
	
	# Check for remaining log files, perhaps from previously completed steps
	my @logs = glob(File::Spec->catfile($self->dir, '*.out.txt'));
	foreach my $log (@logs) {
		push @output, "=== Log file: $log\n";
		if (-z $log) {
			# an empty file
			push @output, "\n";
		}
		else {
			# push log contents to combined output
			my $fh = IO::File->new($log, 'r') or next;
			push @output, <$fh>;
			push @output, "\n";
			$fh->close;
		}
		unlink $log if -e $log;
	}
	
	# print everything out
	my $file = File::Spec->catfile($self->dir, $self->out . "_job_output_logs.txt");
	my $fh = IO::File->new($file, "w");
	foreach (@output) {
		$fh->print($_);
	}
	$fh->close;
	print "\nCombined all job output log files into '$file'\n";
	
	# remove files no longer need
	print "\n\n======= Deleting temporary files\n";
	unless ($self->savebam) {
		foreach my $Job ($self->list_jobs) {
			foreach my $b ($Job->chip_dedup_bams) {
				unlink $b if -e $b;
				$b .= ".bai";
				unlink $b if -e $b;
			}
			foreach my $b ($Job->control_dedup_bams) {
				unlink $b if -e $b;
				$b .= ".bai";
				unlink $b if -e $b;
			}
		}
	}
	if ($self->broad and $self->independent) {
		foreach my $Job ($self->list_jobs) {
			foreach my $f ($Job->rep_gappeaks) {
				my $b = $f;
				$b =~ s/gapped/broad/;
				unlink $b;
			}
		}
	}
	unlink $self->chromofile if $self->chromofile eq 
		File::Spec->catfile($self->dir,"chrom_sizes.temp.txt"); # calculated format
	unlink $self->{progress_file};
}


sub run_organize {
	my $self = shift;
	my $suffix = shift || q();
	return unless ($self->organize);
	return if $self->dryrun;
	print "\n\n======= Moving files into subdirectories\n";
	
	# directories
	my $fragdir  = File::Spec->catfile($self->dir, 'Fragment');
	my $log2dir  = File::Spec->catfile($self->dir, 'Log2FE');
	my $countdir = File::Spec->catfile($self->dir, 'Count');
	my $qdir     = File::Spec->catfile($self->dir, 'QValue');
	my $peakdir  = File::Spec->catfile($self->dir, 'Peaks'. $suffix);
	my $sumitdir = File::Spec->catfile($self->dir, 'PeakSummits' . $suffix);
	my $analdir  = File::Spec->catfile($self->dir, 'Analysis' . $suffix);
	foreach ($fragdir, $log2dir, $countdir, $qdir, $peakdir, $sumitdir, $analdir) {
		make_path($_);
	}
	
	# we're globbing all the files, hope this isn't a problem in case user has 
	# similarly named files in existing directory.
	# Otherwise there's an awful lot of conditional checks for files in every single 
	# ChIPJob object plus general run files....
	
	# log2FE files
	foreach (glob(File::Spec->catfile($self->dir, '*.log2FE.bw')) ) {
		move($_, $log2dir);
	}
	
	# fragment files
	foreach (glob(File::Spec->catfile($self->dir, '*.fragment.bw')) ) {
		move($_, $fragdir);
	}
	foreach (glob(File::Spec->catfile($self->dir, '*.lambda_control.bw')) ) {
		move($_, $fragdir);
	}
	foreach (glob(File::Spec->catfile($self->dir, '*.control_fragment.bw')) ) {
		move($_, $fragdir);
	}
	foreach (glob(File::Spec->catfile($self->dir, '*.fragment.global_mean.bw')) ) {
		move($_, $fragdir);
	}
	
	# qvalue files
	foreach (glob(File::Spec->catfile($self->dir, '*.qvalue.bw')) ) {
		move($_, $qdir);
	}
	
	# count files
	foreach (glob(File::Spec->catfile($self->dir, '*.count.bw')) ) {
		move($_, $countdir);
	}
	
	# peak files
	foreach (glob(File::Spec->catfile($self->dir, '*summit*.bed')) ) {
		move($_, $sumitdir);
	}
	foreach (glob(File::Spec->catfile($self->dir, '*.bed')) ) {
		move($_, $peakdir);
	}
	if ($self->independent) {
		foreach (glob(File::Spec->catfile($self->dir, '*.narrowPeak')) ) {
			move($_, $peakdir);
		}
		if ($self->broad) {
			foreach (glob(File::Spec->catfile($self->dir, '*.gappedPeak')) ) {
				move($_, $peakdir);
			}
		}
	}
	
	# text files
	foreach (glob(File::Spec->catfile($self->dir, '*.txt*')) ) {
		next if $_ =~ /job_output_logs\.txt$/;
		move($_, $analdir);
	}
	
	# independent macs2 analysis files
	if ($self->independent) {
		foreach (glob(File::Spec->catfile($self->dir, '*.xls')) ) {
			move($_, $analdir);
		}
	}
	
	# image files
	if ($self->plot) {
		my $imagedir = File::Spec->catfile($self->dir, 'Plots' . $suffix);
		make_path($imagedir);
		foreach (glob(File::Spec->catfile($self->dir, '*.png')) ) {
			move($_, $imagedir);
		}
	}
	
	# dedup bam files
	if ($self->savebam and $self->dedup) {
		my $bamdir = File::Spec->catfile($self->dir, 'DeDupBam');
		make_path($bamdir);
		foreach (glob(File::Spec->catfile($self->dir, '*.dedup.bam*')) ) {
			move($_, $bamdir);
		}
	}
	
	# saved bedGraph files
	if ($self->savebdg) {
		my $bdgdir = File::Spec->catfile($self->dir, 'BedGraph');
		make_path($bdgdir);
		foreach (glob(File::Spec->catfile($self->dir, '*.bdg')) ) {
			move($_, $bdgdir);
		}
	}
}


