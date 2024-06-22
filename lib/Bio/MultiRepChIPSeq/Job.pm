package Bio::MultiRepChIPSeq::Job;

use strict;
use Carp;
use IO::File;
use File::Spec;
use Data::Dumper;
use List::Util qw(min);
use base 'Bio::MultiRepChIPSeq::options';
use Bio::ToolBox::utility qw(simplify_dataset_name format_with_commas);

our $VERSION = 19.2;

sub new {

	# pass name, comma-list of ChIP files, and comma-list of control files
	my (
		$class, $opts, $job_name, $chip, $control, $chip_scale, $control_scale,
		$chrnorm
	) = @_;

	# bless early to use methods
	my $self = bless {
		opts                => $opts,      # options
		job_name            => $job_name,  # name for this ChIP job
		chip_bams           => [],         # array of initial chip bam file names
		chip_dedup_bams     => [],         # array of deduplicated bam file names
		chip_filter_bams    => [],         # array of filtered bam file names
		chip_use_bams       => [],         # array of final bam file names to use
		chip_count_bw       => [],         # array of ChIP count bigWig file names
		chip_rep_names      => [],         # array of ChIP file base names
		rep_peaks           => [],         # array of ChIP replicate peaks
		rep_gappeaks        => [],         # array of ChIP replicate gapped peaks
		control_bams        => [],         # array of initial control bam file names
		control_dedup_bams  => [],         # array of deduplicated control bam files
		control_filter_bams => [],         # array of filtered control bam files
		control_use_bams    => [],         # array of final control bam files to use
		control_count_bw    => [],         # array of control count bigWig file names
		chip_scale          => [],         # array of scaling factors for chip signal
		control_scale       => [],         # array of scaling factors for control signal
		chr_normfactor      => $chrnorm,   # chromosome-scaling factor
		chip_bw             => undef,      # ChIP signal bigWig file name
		chip_ind_bw         => [],         # individual replicate ChIP signal bigWig names
		chip_bdg            => undef,      # ChIP signal bedGraph file name
		d_control_bdg       => undef,      # d control pileup bedGraph
		d_control_bw        => undef,      # d control pileup bigWig file name
		s_control_bdg       => undef,      # small lambda control bedGraph file name
		l_control_bdg       => undef,      # large lambda control bedGraph file name
		sld_control_file    => undef,      # combined d, small, large lambda file
		lambda_bdg          => undef,      # final lambda control bedGraph file name
		lambda_bw           => undef,      # final lambda control bigWig file name
		fe_bdg              => undef,      # fold enrichment bedGraph file name
		logfe_bw            => undef,      # log2 fold enrichment bigWig file name
		qvalue_bdg          => undef,      # q-value bedGraph file name
		qvalue_bw           => undef,      # q-value bigWig file name
		repmean_peak        => undef,      # mean-replicate peak file name
		repmean_gappeak     => undef,      # mean-replicate gapped peak file name
		repmerge_peak       => undef,      # merged-replicate peak file name
		repmerge_gappeak    => undef,      # merged-replicate gapped peak file name
	}, $class;

	# base namepath
	my $namepath = File::Spec->catfile( $self->dir, $job_name );

	## check the ChIP files
	if ( $chip =~ /\. (?: bw | bigwig )$/xi ) {
		$self->chip_bw($chip);
		$self->crash("only one ChIP bigWig file is allowed per experiment!\n")
			if $chip =~ /,/;
		if ( $self->independent ) {
			$self->crash("cannot use independent flag with bigwig files!\n");
		}
		$self->chip_bdg("$namepath.fragment.bdg");
		$self->qvalue_bdg("$namepath.qvalue.bdg");
		$self->qvalue_bw("$namepath.qvalue.bw");
		$self->fe_bdg("$namepath.FE.bdg");
		$self->logfe_bw("$namepath.log2FE.bw");
		$self->repmean_peak("$namepath.narrowPeak");
		$self->repmean_gappeak("$namepath.gappedPeak");
	}
	elsif ( $chip =~ /\.bam$/i ) {
		my @bams = split( /,/, $chip );
		$self->chip_bams(@bams);
		$self->chip_bdg("$namepath.fragment.bdg");
		$self->chip_bw("$namepath.fragment.bw");
		$self->qvalue_bdg("$namepath.qvalue.bdg");
		$self->qvalue_bw("$namepath.qvalue.bw");
		$self->fe_bdg("$namepath.FE.bdg");
		$self->logfe_bw("$namepath.log2FE.bw");

		if ( $self->independent ) {
			$self->repmean_peak("$namepath.rep_mean.narrowPeak");
			$self->repmean_gappeak("$namepath.rep_mean.gappedPeak");
			$self->repmerge_peak("$namepath.rep_merge.bed");
			$self->repmerge_gappeak("$namepath.rep_merge.broad.bed");
		}
		else {
			$self->repmean_peak("$namepath.narrowPeak");
			$self->repmean_gappeak("$namepath.gappedPeak");
		}
		if ($chip_scale) {
			$self->chip_scale( split( /,/, $chip_scale ) );
			$self->crash("unequal scale factors and bam files!\n")
				if scalar( $self->chip_scale ) != scalar(@bams);
		}

		# generate count bw and dedup bam file names
		foreach my $bam (@bams) {
			my ( undef, undef, $fname ) = File::Spec->splitpath($bam);
			$fname = simplify_dataset_name($fname);    # strip extension and cleanup
			$self->chip_rep_names($fname);
			my $base = File::Spec->catfile( $self->dir, $fname );

			# count bigWig
			$self->chip_count_bw("$base.count.bw");

			# dedup and filter bams
			$self->chip_dedup_bams("$base.dedup.bam");
			$self->chip_filter_bams("$base.filter.bam");

			# if running independent peak calls
			if ( $self->independent ) {

				# add peak and coverage files
				$self->rep_peaks( $base . '_peaks.narrowPeak' );
				$self->rep_gappeaks( $base . '_peaks.gappedPeak' );
				$self->chip_ind_bw("$base.fragment.bw");
			}
		}
	}

	## check the control files
	if ( $control =~ /^Custom \- Universal \- (.+)$/x ) {

		# this is just a ChIP Job using a common universal reference file
		# the actual control name is embedded in the name, so extract it
		# we only need to know the bedgraph file name here
		$self->lambda_bdg( File::Spec->catfile( $opts->{dir}, $1 . '.bdg' ) );
	}
	elsif ( $control =~ /\.(?: bw | bigwig )$/xi ) {

		# pre-generated reference bigWig file
		$self->crash("only one control bigWig file is allowed per experiment!\n")
			if $chip =~ /,/;
		$self->lambda_bw($control);
		$self->lambda_bdg( $namepath . '.bdg' );
	}
	elsif ( $control =~ /\.bam$/i ) {

		# we have control bams to process
		my @bams = split( /,/, $control );
		$self->control_bams(@bams);
		if ( $self->lambda ) {
			my $control_base = $namepath;
			if ( $namepath =~ /\.lambda_control/x ) {

				# the lambda_control extension is already appended to the name
				# as is the case with universal controls
				$self->lambda_bdg( $namepath . '.bdg' );
				$self->lambda_bw( $namepath . '.bw' );

				# strip the lambda_control bit for the other intermediate files
				$control_base =~ s/\.lambda_control//x;
			}
			else {
				# append the lambda_control extension
				$self->lambda_bdg( $namepath . '.lambda_control.bdg' );
				$self->lambda_bw( $namepath . '.lambda_control.bw' );
			}
			$self->d_control_bdg("$control_base.dlocal.bdg");
			$self->d_control_bw("$control_base.control_fragment.bw");
			$self->s_control_bdg("$control_base.slocal.bdg") if $self->slocal;
			$self->l_control_bdg("$control_base.llocal.bdg") if $self->llocal;
			$self->sld_control_file("$control_base.sldlocal.txt");
		}
		else {
			# even with no lambda, we recycle the hash entries for simplicity
			if ( $namepath =~ /\.control_fragment/x ) {

				# the fragment control extension is already appended to the name
				# as is the case with universal controls
				$self->lambda_bdg( $namepath . '.bdg' );
				$self->lambda_bw( $namepath . '.bw' );
			}
			else {
				# append the fragment control extension
				$self->lambda_bdg( $namepath . '.control_fragment.bdg' );
				$self->lambda_bw( $namepath . '.control_fragment.bw' );
			}
		}

		if ($control_scale) {
			$self->control_scale( split( /,/, $control_scale ) );
			$self->crash("unequal scale factors and bam files!\n")
				if scalar(@bams) != scalar( $self->control_scale );
		}

		# generate count bw  and dedup file names
		foreach my $bam (@bams) {
			my ( undef, undef, $fname ) = File::Spec->splitpath($bam);
			$fname =~ s/\.bam$//i;    # strip extension
			my $base = File::Spec->catfile( $self->dir, $fname );

			# count bigWig
			$self->control_count_bw("$base.count.bw");

			# dedup and filter bams
			$self->control_dedup_bams("$base.dedup.bam");
			$self->control_filter_bams("$base.filter.bam");
		}
	}
	else {
		# must be just a chip without corresponding control
		$self->lambda_bdg("$namepath.fragment.global_mean.bdg");
		$self->lambda_bw("$namepath.fragment.global_mean.bw");
	}

	return $self;
}

sub crash {
	my ( $self, $message ) = @_;
	printf STDERR "ChIP Job object:\n%s\n", Dumper($self);
	confess $message;
}

sub job_name {
	my $self = shift;
	$self->{job_name} = shift if @_;
	return $self->{job_name};
}

sub chip_bams {
	my $self = shift;
	push @{ $self->{chip_bams} }, @_ if @_;
	return @{ $self->{chip_bams} };
}

sub control_bams {
	my $self = shift;
	push @{ $self->{control_bams} }, @_ if @_;
	return @{ $self->{control_bams} };
}

sub chip_dedup_bams {
	my $self = shift;
	push @{ $self->{chip_dedup_bams} }, @_ if @_;
	return @{ $self->{chip_dedup_bams} };
}

sub chip_filter_bams {
	my $self = shift;
	push @{ $self->{chip_filter_bams} }, @_ if @_;
	return @{ $self->{chip_filter_bams} };
}

sub chip_use_bams {
	my $self = shift;
	push @{ $self->{chip_use_bams} }, @_ if @_;
	return @{ $self->{chip_use_bams} };
}

sub control_dedup_bams {
	my $self = shift;
	push @{ $self->{control_dedup_bams} }, @_ if @_;
	return @{ $self->{control_dedup_bams} };
}

sub control_filter_bams {
	my $self = shift;
	push @{ $self->{control_filter_bams} }, @_ if @_;
	return @{ $self->{control_filter_bams} };
}

sub control_use_bams {
	my $self = shift;
	push @{ $self->{control_use_bams} }, @_ if @_;
	return @{ $self->{control_use_bams} };
}

sub chip_scale {
	my $self = shift;
	push @{ $self->{chip_scale} }, @_ if @_;
	return @{ $self->{chip_scale} };
}

sub control_scale {
	my $self = shift;
	push @{ $self->{control_scale} }, @_ if @_;
	return @{ $self->{control_scale} };
}

sub chr_normfactor {
	my $self = shift;
	$self->{chr_normfactor} = shift if @_;
	return $self->{chr_normfactor};
}

sub chip_rep_names {
	my $self = shift;
	push @{ $self->{chip_rep_names} }, @_ if @_;
	return @{ $self->{chip_rep_names} };
}

sub rep_peaks {
	my $self = shift;
	push @{ $self->{rep_peaks} }, @_ if @_;
	return @{ $self->{rep_peaks} };
}

sub rep_gappeaks {
	my $self = shift;
	push @{ $self->{rep_gappeaks} }, @_ if @_;
	return @{ $self->{rep_gappeaks} };
}

sub chip_bw {
	my $self = shift;
	$self->{chip_bw} = shift if @_;
	return $self->{chip_bw};
}

sub chip_bdg {
	my $self = shift;
	$self->{chip_bdg} = shift if @_;
	return $self->{chip_bdg};
}

sub chip_ind_bw {
	my $self = shift;
	push @{ $self->{chip_ind_bw} }, @_ if @_;
	return @{ $self->{chip_ind_bw} };
}

sub chip_count_bw {
	my $self = shift;
	push @{ $self->{chip_count_bw} }, @_ if @_;
	return @{ $self->{chip_count_bw} };
}

sub control_count_bw {
	my $self = shift;
	push @{ $self->{control_count_bw} }, @_ if @_;
	return @{ $self->{control_count_bw} };
}

sub d_control_bdg {
	my $self = shift;
	$self->{d_control_bdg} = shift if @_;
	return $self->{d_control_bdg};
}

sub d_control_bw {
	my $self = shift;
	$self->{d_control_bw} = shift if @_;
	return $self->{d_control_bw};
}

sub s_control_bdg {
	my $self = shift;
	$self->{s_control_bdg} = shift if @_;
	return $self->{s_control_bdg};
}

sub l_control_bdg {
	my $self = shift;
	$self->{l_control_bdg} = shift if @_;
	return $self->{l_control_bdg};
}

sub sld_control_file {
	my $self = shift;
	$self->{sld_control_file} = shift if @_;
	return $self->{sld_control_file};
}

sub lambda_bdg {
	my $self = shift;
	$self->{lambda_bdg} = shift if @_;
	return $self->{lambda_bdg};
}

sub lambda_bw {
	my $self = shift;
	$self->{lambda_bw} = shift if @_;
	return $self->{lambda_bw};
}

sub fe_bdg {
	my $self = shift;
	$self->{fe_bdg} = shift if @_;
	return $self->{fe_bdg};
}

sub logfe_bw {
	my $self = shift;
	$self->{logfe_bw} = shift if @_;
	return $self->{logfe_bw};
}

sub qvalue_bdg {
	my $self = shift;
	$self->{qvalue_bdg} = shift if @_;
	return $self->{qvalue_bdg};
}

sub qvalue_bw {
	my $self = shift;
	$self->{qvalue_bw} = shift if @_;
	return $self->{qvalue_bw};
}

sub repmean_peak {
	my $self = shift;
	$self->{repmean_peak} = shift if @_;
	return $self->{repmean_peak};
}

sub repmean_gappeak {
	my $self = shift;
	$self->{repmean_gappeak} = shift if @_;
	return $self->{repmean_gappeak};
}

sub repmerge_peak {
	my $self = shift;
	$self->{repmerge_peak} = shift if @_;
	return $self->{repmerge_peak};
}

sub repmerge_gappeak {
	my $self = shift;
	$self->{repmerge_gappeak} = shift if @_;
	return $self->{repmerge_gappeak};
}

sub generate_dedup_commands {
	my $self      = shift;
	my $name2done = shift;
	unless ( $self->bamdedup_app =~ /\w+/ or $self->dryrun ) {
		croak "no bam_partial_dedup.pl application in path!\n";
	}

	# first collect all the bam files
	my @bamfiles;
	if ( $self->chip_bams ) {
		my @b  = $self->chip_bams;
		my @db = $self->chip_dedup_bams;
		for my $i ( 0 .. $#b ) {
			push @bamfiles, [
				$i,
				'chip',
				$b[$i],    # input bam
				$db[$i]    # output dedup bam
			];
		}
	}
	if ( $self->control_bams ) {
		my @b  = $self->control_bams;
		my @db = $self->control_dedup_bams;
		for my $i ( 0 .. $#b ) {
			push @bamfiles, [
				$i,
				'control',
				$b[$i],    # input bam
				$db[$i]    # output dedup bam
			];
		}
	}

	# generate the commands
	my @commands;
	foreach my $set (@bamfiles) {
		my $in  = $set->[2];
		my $out = $set->[3];

		# check if this has been done, should only matter with shared control files
		if ( exists $name2done->{$in} ) {

			# this file has already been done, but we need to update the name
			# this is quite ugly, but I want to be able to do an exact replacement
			if ( $set->[1] eq 'control' ) {
				$self->{control_dedup_bams}->[ $set->[0] ] = $name2done->{$in};
			}
			elsif ( $set->[1] eq 'chip' ) {

				# this is likely a mistake or duplicate
				$self->{chip_dedup_bams}->[ $set->[0] ] = $name2done->{$in};
			}
			next;
		}

		# generate command
		my $command = sprintf "%s --in %s --out %s --cpu %s --qual %d ",
			$self->bamdedup_app || 'bam_partial_dedup.pl',
			$in,
			$out,
			$self->cpu,
			$self->mapq;
		if ( $self->paired or $self->deduppair ) {
			$command .= "--pe ";
		}
		if ( $self->maxdepth and $self->maxdepth == 1 ) {

			# no other deduplication options need to be set
			$command .= "--max 1 ";
		}
		else {
			# set random subsampling, maximum duplicates, and/or optical
			if ( $self->dupfrac > 0 ) {
				$command .= sprintf "--seed 1 --frac %s ", $self->dupfrac;
			}
			if ( $self->maxdepth and $self->maxdepth > 1 ) {
				$command .= sprintf "--max %s ", $self->maxdepth;
			}
			if ( $self->optdist ) {
				$command .= sprintf "--optical --distance %s ", $self->optdist;
			}
		}
		if ( $self->blacklist and $self->blacklist ne 'none' ) {
			$command .= sprintf "--blacklist %s ", $self->blacklist;
		}
		if ( $self->chrskip ) {
			$command .= sprintf "--chrskip \'%s\' ", $self->chrskip;
		}
		my $log = $out;
		$log =~ s/\.bam$/.out.txt/i;
		$command .= " 2>&1 > $log";
		push @commands, [ $command, $out, $log ];
		$name2done->{$in} = $out;    # remember that this has been done
	}

	return @commands;
}

sub generate_bam_filter_commands {
	my $self        = shift;
	my $name2done   = shift;
	my $filter_file = shift;
	unless ( $self->bedtools_app =~ /\w+/ or $self->dryrun ) {
		croak "no bedtools application in path!\n";
	}
	unless ( $self->samtools_app =~ /\w+/ or $self->dryrun ) {
		croak "no samtools application in path!\n";
	}

	# first collect all the bam files
	my @bamfiles;
	if ( $self->chip_bams ) {
		my @b  = $self->chip_bams;
		my @fb = $self->chip_filter_bams;
		for my $i ( 0 .. $#b ) {
			push @bamfiles, [
				$i,
				'chip',
				$b[$i],    # input bam
				$fb[$i]    # output filter bam
			];
		}
	}
	if ( $self->control_bams ) {
		my @b  = $self->control_bams;
		my @db = $self->control_filter_bams;
		for my $i ( 0 .. $#b ) {
			push @bamfiles, [
				$i,
				'control',
				$b[$i],    # input bam
				$db[$i]    # output filter bam
			];
		}
	}

	# generate samtool filter options
	my $filter;
	if ( $self->paired or $self->deduppair ) {

		# exclude filter: UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
		# include filter: PROPER_PAIR
		$filter = "-f 2 -F 3852 ";
	}
	else {
		# exclude filter: UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
		$filter = "-F 3844 ";
	}
	if ( $self->mapq ) {
		$filter .= sprintf "-q %d ", $self->mapq;
	}

	# generate the commands
	my @commands;
	foreach my $set (@bamfiles) {
		my $in  = $set->[2];
		my $out = $set->[3];

		# check if this has been done, should only matter with shared control files
		if ( exists $name2done->{$in} ) {

			# this file has already been done, but we need to update the name
			# this is quite ugly, but I want to be able to do an exact replacement
			if ( $set->[1] eq 'control' ) {
				$self->{control_filter_bams}->[ $set->[0] ] = $name2done->{$in};
			}
			elsif ( $set->[1] eq 'chip' ) {

				# this is likely a mistake or duplicate
				$self->{chip_filter_bams}->[ $set->[0] ] = $name2done->{$in};
			}
			next;
		}

		# generate command
		my $sam_command = sprintf "%s view $filter --threads %d --bam --output %s",
			$self->samtools_app,
			$self->cpu,
			$out;
		my $command;
		if ($filter_file) {

			# bedtools piped to samtools
			$command = sprintf "%s intersect -v -ubam -a %s -b %s | %s - ",
				$self->bedtools_app,
				$in,
				$filter_file,
				$sam_command;
		}
		else {
			# samtools alone
			$command = sprintf "%s %s ", $sam_command, $in;
		}
		my $log = $out;
		$log =~ s/\.bam$/.out.txt/i;
		$command .= "2>&1 > $log";
		push @commands, [ $command, $out, $log ];
		$name2done->{$in} = $out;    # remember that this has been done
	}

	return @commands;
}

sub find_new_bams {
	my $self = shift;

	# ChIP bams
	if ( $self->chip_bams ) {
		printf " Checking ChIP bams for %s:\n", $self->job_name;
		my @b  = $self->chip_bams;
		my @db = $self->chip_dedup_bams;
		my @fb = $self->chip_filter_bams;
		for my $i ( 0 .. $#b ) {
			my $in   = $b[$i];
			my $out1 = $db[$i];
			my $out2 = $fb[$i];
			if ( -e $out1 and -s _ ) {
				$self->chip_use_bams($out1);
				printf "  Found $out1\n";
			}
			elsif ( -e $out2 and -s _ ) {
				$self->chip_use_bams($out2);
				printf "  Found $out2\n";
			}
			else {
				$self->chip_use_bams($in);
				printf "  Using $in\n";
			}
		}
	}

	# Control bams
	if ( $self->control_bams ) {
		printf " Checking control bams for %s:\n", $self->job_name;
		my @b  = $self->control_bams;
		my @db = $self->control_dedup_bams;
		my @fb = $self->control_filter_bams;
		for my $i ( 0 .. $#b ) {
			my $in   = $b[$i];
			my $out1 = $db[$i];
			my $out2 = $fb[$i];
			if ( -e $out1 and -s _ ) {
				$self->control_use_bams($out1);
				printf "  Found $out1\n";
			}
			elsif ( -e $out2 and -s _ ) {
				$self->control_use_bams($out2);
				printf "  Found $out2\n";
			}
			else {
				$self->control_use_bams($in);
				printf "  Using $in\n";
			}
		}
	}
}

sub generate_bam2wig_frag_commands {
	my $self      = shift;
	my $name2done = shift;
	unless ( $self->bam2wig_app =~ /\w+/ or $self->dryrun ) {
		croak "no bam2wig.pl application in path!\n";
	}
	my @commands;

	### ChIP bams
	if ( $self->chip_use_bams ) {

		# we have bam files to convert to bw

		# base fragment command
		my $frag_base = sprintf
			"--qual %s --nosecondary --noduplicate --nosupplementary --bin %s --cpu %s ",
			$self->mapq, $self->chipbin, $self->cpu;

		# paired or single options
		if ( $self->paired ) {
			$frag_base .= sprintf "--pe --span --minsize %s --maxsize %s ",
				$self->minsize, $self->maxsize;
		}
		else {
			$frag_base .= sprintf "--extend --extval %s ", $self->fragsize;
			if ( $self->shiftsize ) {
				$frag_base .= sprintf "--shiftval %s ", $self->shiftsize;
			}
		}

		# additional filtering
		if ( $self->fraction ) {
			$frag_base .= "--fraction ";
		}
		if ( $self->blacklist and $self->blacklist ne 'none' ) {
			$frag_base .= sprintf "--blacklist %s ", $self->blacklist;
		}
		if ( $self->chrskip ) {
			$frag_base .= sprintf "--chrskip \'%s\' ", $self->chrskip;
		}

		# scaling
		if ( $self->chip_scale ) {
			$frag_base .= sprintf "--scale %s ", join( ',', $self->chip_scale );
		}
		else {
			# standard scaling
			$frag_base .= "--rpm ";
		}

		# chromosome-specific scaling
		if ( $self->chr_normfactor ) {
			$frag_base .= sprintf "--chrnorm %s --chrapply %s ",
				$self->chr_normfactor, $self->chrapply;
		}

		# finish fragment command
		my $frag_command = sprintf
			"%s --out %s %s --mean --bdg --in %s",
			$self->bam2wig_app || 'bam2wig.pl',
			$self->chip_bdg,
			$frag_base,
			join( ',', $self->chip_use_bams );
		my $log = $self->chip_bdg;
		$log =~ s/bdg$/bam2wig.out.txt/;
		$frag_command .= " 2>&1 > $log";
		push @commands, [ $frag_command, $self->chip_bdg, $log ];

		### Independent ChIP bams
		if ( $self->independent ) {
			my @bams = $self->chip_use_bams;
			for my $i ( 0 .. $#bams ) {
				my $out     = ( $self->chip_ind_bw )[$i];
				my $command = sprintf
					"%s --out %s %s --bw --bwapp %s --in %s",
					$self->bam2wig_app || 'bam2wig.pl',
					$out,
					$frag_base,
					$self->wig2bw_app || 'wigToBigWig',
					$bams[$i];
				$log = $out;
				$log =~ s/bw/bam2wig.out.txt/;
				$command .= " 2>&1 > $log";
				push @commands, [ $command, $out, $log ];
			}
		}
	}

	### Control bams
	my $control_bam_string = join( ',', $self->control_use_bams );

	# first check if we've processed this control yet
	if ( exists $name2done->{$control_bam_string} ) {

		# we have, so change the lambda bedgraph name to the file that's been used
		# this will look funny, because it will look like we're using another chip's
		# control file, but that's the way it will be
		$self->lambda_bdg( $name2done->{$control_bam_string} );
	}

	# process control bams into a chromatin bias lambda-control track
	elsif ( scalar( $self->control_use_bams )
		and $self->lambda
		and not exists $name2done->{$control_bam_string} )
	{

		# base fragment command
		my $frag_command = sprintf
"--qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --mean --bdg ",
			$self->mapq, $self->cpu;

		# general paired options, restrict size for all
		if ( $self->paired ) {
			$frag_command .= sprintf "--pe --minsize %s --maxsize %s ",
				$self->minsize, $self->maxsize;
		}

		# additional filters
		if ( $self->fraction ) {
			$frag_command .= "--fraction ";
		}
		if ( $self->blacklist and $self->blacklist ne 'none' ) {
			$frag_command .= sprintf "--blacklist %s ", $self->blacklist;
		}
		if ( $self->chrskip ) {
			$frag_command .= sprintf "--chrskip \'%s\' ", $self->chrskip;
		}

		# chromosome-specific scaling
		if ( $self->chr_normfactor ) {
			$frag_command .= sprintf "--chrnorm %s --chrapply %s ",
				$self->chr_normfactor, $self->chrapply;
		}

		## now duplicate the base command for each size lambda control

		# d control, use extend or paired span just like ChIP
		my $command1 = $frag_command;
		if ( $self->paired ) {
			$command1 .= "--span ";
		}
		else {
		 # treat d just like ChIP, including extend and shift, this may deviate from Macs2
			$command1 .= sprintf "--extend --extval %s --shiftval %0.0f ",
				$self->fragsize, $self->shiftsize;
		}
		if ( $self->control_scale ) {
			$command1 .= sprintf "--scale %s ", join( ',', $self->control_scale );
		}
		else {
			$command1 .= "--rpm ";
		}

		# regenerate command 1 with program and input/output files
		$command1 = sprintf
			"%s --out %s %s --bin %s --in %s ",
			$self->bam2wig_app || 'bam2wig.pl',
			$self->d_control_bdg,
			$command1,
			$self->chipbin,
			$control_bam_string;
		my $log = $self->d_control_bdg;
		$log =~ s/bdg$/bam2wig.out.txt/;
		$command1 .= " 2>&1 > $log";
		push @commands, [ $command1, $self->d_control_bdg, $log ];

		# small local lambda, extend both directions, scaled to compensate for length
		if ( $self->slocal ) {
			my $command2 = $frag_command;
			my $scale    = sprintf "%.4f", $self->fragsize / $self->slocal;
			if ( $self->control_scale ) {

				# user provided scale, multiply this with the lambda size scale
				my @scales =
					map { sprintf "%.4f", $_ * $scale } ( $self->control_scale );
				$scale = join( ',', @scales );    # replace
			}
			else {
				# standard scaling
				$command2 .= "--rpm ";
			}

			# regenerate command 2 with program and input/output files
			$command2 = sprintf
				"%s --out %s %s --cspan --extval %s --scale %s --bin %s --in %s ",
				$self->bam2wig_app || 'bam2wig.pl',
				$self->s_control_bdg,
				$command2,
				$self->slocal,
				$scale,
				$self->slocalbin,
				$control_bam_string;
			$log = $self->s_control_bdg;
			$log =~ s/bdg$/bam2wig.out.txt/;
			$command2 .= " 2>&1 > $log";
			push @commands, [ $command2, $self->s_control_bdg, $log ];
		}

		# large local lambda, extend both directions, scaled to compensate for length
		if ( $self->llocal ) {
			my $command3 = $frag_command;
			my $scale    = sprintf "%.4f", $self->fragsize / $self->llocal;
			if ( $self->control_scale ) {

				# user provided scale, multiply this with the lambda size scale
				my @scales =
					map { sprintf "%.4f", $_ * $scale } ( $self->control_scale );
				$scale = join( ',', @scales );    # replace
			}
			else {
				# standard scaling
				$command3 .= "--rpm ";
			}

			# regenerate command 3 with program and input/output files
			$command3 = sprintf
				"%s --out %s %s --cspan --extval %s --scale %s --bin %s --in %s ",
				$self->bam2wig_app || 'bam2wig.pl',
				$self->l_control_bdg,
				$command3,
				$self->llocal,
				$scale,
				$self->llocalbin,
				$control_bam_string;
			$log = $self->l_control_bdg;
			$log =~ s/bdg$/bam2wig.out.txt/;
			$command3 .= " 2>&1 > $log";
			push @commands, [ $command3, $self->l_control_bdg, $log ];
		}

	  # record that we've done this bam
	  # we store the lambda bedgraph name because it might be reused again for another job
		$name2done->{$control_bam_string} = $self->lambda_bdg;
	}

	# skipping chromatin-bias lambda control track, use control track as is
	elsif ( scalar( $self->control_use_bams )
		and not $self->lambda
		and not exists $name2done->{$control_bam_string} )
	{

		# fragment command
		my $frag_command = sprintf
"%s --out %s --qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --mean --bdg ",
			$self->bam2wig_app || 'bam2wig.pl',
			$self->lambda_bdg,
			$self->mapq,
			$self->cpu;

		# single or paired options
		if ( $self->paired ) {
			$frag_command .= sprintf "--span --pe --minsize %s --maxsize %s --bin %s ",
				$self->minsize,
				$self->maxsize,
				$self->chipbin;
		}
		else {
			$frag_command .= sprintf "--extend --extval %s --shiftval %0.0f --bin %s ",
				$self->fragsize,
				$self->shiftsize,
				$self->chipbin;
		}

		# scaling for the fragment command only
		if ( $self->control_scale ) {
			$frag_command .= sprintf "--scale %s ", join( ',', ( $self->control_scale ) );
		}
		else {
			# standard scaling
			$frag_command .= "--rpm ";
		}

		# additional filters
		if ( $self->blacklist and $self->blacklist ne 'none' ) {
			$frag_command .= sprintf "--blacklist %s ", $self->blacklist;
		}
		if ( $self->chrskip ) {
			$frag_command .= sprintf "--chrskip \'%s\' ", $self->chrskip;
		}

		# chromosome-specific scaling
		if ( $self->chr_normfactor ) {
			$frag_command .= sprintf "--chrnorm %s --chrapply %s ",
				$self->chr_normfactor, $self->chrapply;
		}

		# finish command
		$frag_command .= "--in $control_bam_string ";
		my $log = $self->lambda_bdg;
		$log =~ s/bdg$/bam2wig.out.txt/;
		$frag_command .= " 2>&1 > $log";
		push @commands, [ $frag_command, $self->lambda_bdg, $log ];

	  # record that we've done this bam
	  # we store the lambda bedgraph name because it might be reused again for another job
		$name2done->{$control_bam_string} = $self->lambda_bdg;
	}

	# finished
	return @commands;
}

sub generate_bam2wig_count_commands {
	my $self      = shift;
	my $name2done = shift;
	unless ( $self->bam2wig_app =~ /\w+/ or $self->dryrun ) {
		croak "no bam2wig.pl application in path!\n";
	}
	unless ( $self->wig2bw_app =~ /\w+/ or $self->dryrun ) {
		croak "no wigToBigWig application in path!\n";
	}
	my @commands;

	### ChIP bams
	if ( $self->chip_use_bams ) {

		# we have bam files to count

		# base count command
		my $count_command = sprintf
"--qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --bw --bwapp %s ",
			$self->mapq,
			$self->cpu,
			$self->wig2bw_app || 'wigToBigWig';

		# paired or single options
		if ( $self->paired ) {
			$count_command .= sprintf "--pe --mid --minsize %s --maxsize %s ",
				$self->minsize, $self->maxsize;
		}
		elsif ( $self->atac ) {
			# Special options for ATAC cutsite point data.
			# Recommended options are to shift by unequal amounts based on strand
			# +5/-4 for F and R reads, but bam2wig isn't set up to do handle that
			# so we simply shift by 5 bp (strand inherently handled).
			# Also bam2wig has special ends mode to make this a little more efficient
			# but doesn't currently handle shift, so just simply run as single end.
			# The bam file should be filtered anyway to remove straggling singletons.
			$count_command .= sprintf "--start --shiftval 5 ";
		}
		else {
			# default single-end
			$count_command .= sprintf "--start --shiftval %0.0f ",
				( $self->fragsize / 2 ) + $self->shiftsize;
		}

		# additional filtering
		if ( $self->fraction ) {
			$count_command .= "--fraction ";
		}
		if ( $self->blacklist and $self->blacklist ne 'none' ) {
			$count_command .= sprintf "--blacklist %s ", $self->blacklist;
		}
		if ( $self->chrskip ) {
			$count_command .= sprintf "--chrskip \'%s\' ", $self->chrskip;
		}

		# finish count commands
		my @bams = $self->chip_use_bams;
		for my $i ( 0 .. $#bams ) {
			my $command = $count_command;

			# add scaling as necessary
			unless ( $self->rawcounts ) {
				if ( $self->chip_scale ) {
					$command .= sprintf "--scale %s ", ( $self->chip_scale )[$i];
				}
				else {
					# we scale the count back up to target depth
					$command .= sprintf "--rpm --scale %s ", $self->targetdepth;
				}
				if ( $self->chr_normfactor ) {

					# chromosome-specific scaling
					$count_command .= sprintf "--chrnorm %s --chrapply %s ",
						$self->chr_normfactor, $self->chrapply;
				}
				$count_command .= "--format 3 ";    # always limit digits
			}

			# regenerate command with file names
			my $out = ( $self->chip_count_bw )[$i];
			$command = sprintf
				"%s --out %s %s --in %s ",
				$self->bam2wig_app || 'bam2wig.pl',
				$out,
				$command,
				$bams[$i];

			# output log
			my $log = $out;
			$log =~ s/bw$/bam2wig.out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [ $command, $out, $log ];
		}
	}

	### Control bams
	# we only process this if we haven't done them yet
	my $control_bam_string = join( ',', ( $self->control_use_bams ) );
	if ( $control_bam_string and not exists $name2done->{$control_bam_string} ) {

		# base count command
		my $count_command = sprintf
"--qual %s --nosecondary --noduplicate --nosupplementary --cpu %s --bw --bwapp %s ",
			$self->mapq,
			$self->cpu,
			$self->wig2bw_app || 'wigToBigWig';

		# general paired options, restrict size for all
		if ( $self->paired ) {
			$count_command .= sprintf "--pe --mid --minsize %s --maxsize %s ",
				$self->minsize, $self->maxsize;
		}
		else {
			$count_command .= sprintf "--start --shiftval %0.0f ",
				( $self->fragsize / 2 ) + $self->shiftsize;
		}

		# additional filters
		if ( $self->fraction ) {
			$count_command .= "--fraction ";
		}
		if ( $self->blacklist and $self->blacklist ne 'none' ) {
			$count_command .= sprintf "--blacklist %s ", $self->blacklist;
		}
		if ( $self->chrskip ) {
			$count_command .= sprintf "--chrskip \'%s\' ", $self->chrskip;
		}

		# add count commands
		my @bams = $self->control_use_bams;
		for my $i ( 0 .. $#bams ) {
			my $command = $count_command;

			# add scaling as necessary
			unless ( $self->rawcounts ) {
				if ( $self->control_scale ) {
					$command .= sprintf "--scale %s ", ( $self->control_scale )[$i];
				}
				else {
					# we scale the count back up to target depth
					$command .= sprintf "--rpm --scale %s ", $self->targetdepth;
				}
				if ( $self->chr_normfactor ) {

					# chromosome-specific scaling
					$count_command .= sprintf "--chrnorm %s --chrapply %s ",
						$self->chr_normfactor, $self->chrapply;
				}
				$command .= "--format 3 ";
			}

			# regenerate command with file names
			my $out = ( $self->control_count_bw )[$i];
			$command = sprintf "%s --out %s %s --in %s ",
				$self->bam2wig_app || 'bam2wig.pl',
				$out,
				$command,
				$bams[$i];

			# output log
			my $log = $out;
			$log =~ s/bw$/bam2wig.out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [ $command, $out, $log ];
		}

		# record that we've done this bam
		$name2done->{$control_bam_string} = 1;
	}

	# finished
	return @commands;
}

sub generate_lambda_control_commands {
	my $self      = shift;
	my $name2done = shift;
	return if exists $name2done->{ $self->lambda_bdg };    # already done

	# check whether no reference control was provided at all
	if (    not scalar( $self->control_use_bams )
		and $self->lambda_bdg =~ /global_mean \.bdg $/x
		and not -e $self->lambda_bdg )
	{
	# no control bam files at all, need to use expected mean
	# generate this from the existing ChIP fragment bigWig file
	# no output file needed to be specified, it will automatically append .global_mean.bdg
		unless ( $self->meanbdg_app =~ /\w+/ or $self->dryrun ) {
			croak "no generate_mean_bedGraph.pl application in path!\n";
		}
		my $log = $self->lambda_bdg;
		$log =~ s/bdg/meanbdg.out.txt/;
		my $chipfile;    # we can use either bw or bdg
						 # but prefer bdg to match the original chromosome order
						 # without having to resort to separate expensive reordering
						 # since our bw adapter dumps chromosomes in random order - sigh
						 # this should in most cases default to the bdg

		if ( -e $self->chip_bdg ) {
			$chipfile = $self->chip_bdg;
		}
		elsif ( -e $self->chip_bw ) {
			$chipfile = $self->chip_bw;
		}
		my $command = sprintf
			"%s --in %s --out %s --genome %d 2>&1 > $log",
			$self->meanbdg_app || 'generate_mean_bedGraph.pl',
			$chipfile,
			$self->lambda_bdg,
			$self->genome || 0;    # value of zero will default to empirical
		$name2done->{ $self->lambda_bdg } = 1;
		return [ $command, $self->lambda_bdg, $log ];
	}

	# Proceed with generating lambda bedGraph
	unless ( $self->macs_app =~ /\w+/ or $self->dryrun ) {
		croak "no MACS2 application in path!\n";
	}
	unless ( $self->bedtools_app =~ /\w+/ or $self->dryrun ) {
		croak "no bedtools application in path!\n";
	}
	unless ( $self->data2wig_app =~ /\w+/ or $self->dryrun ) {
		croak "no data2wig.pl application in path!\n";
	}
	my $dfile = $self->d_control_bdg;
	my $sfile = $self->s_control_bdg;
	my $lfile = $self->l_control_bdg;
	return unless ($dfile);    # controls with lambda will always have  d_control_bdg
							   # ChIP jobs and non-lambda controls will not
	unless ( $self->dryrun ) {
		$self->crash("no d control bedGraph file '$dfile'!\n") if ( not -e $dfile );
		$self->crash("no small control bedGraph file '$sfile'!\n")
			if ( $sfile and not -e $sfile );
		$self->crash("no large control bedGraph file '$lfile'!\n")
			if ( $lfile and not -e $lfile );
	}

	# generate background lambda bedgraph
	# go ahead and do this immediately because it's so simple
	my $background_bdg;
	unless ( $self->dryrun ) {
		my $background = sprintf "%.4f", ( 1_000_000 * $self->fragsize ) / $self->genome;
		printf " Calculating background for %s as $background\n", $self->job_name;
		$background_bdg = $self->lambda_bdg;
		$background_bdg =~ s/lambda_control/background/;
		my $chromofile = $self->chromofile;
		my $infh       = IO::File->new( $chromofile, 'r' )
			or die "unable to open chromosome file '$chromofile'!\n";
		my $outfh = IO::File->new( $background_bdg, "w" );

		while ( my $line = $infh->getline ) {
			chomp $line;
			my ( $chr, $end ) = split /\s/, $line;
			$outfh->printf( "%s\t0\t%s\t%s\n", $chr, $end, $background );
		}
		$infh->close;
		$outfh->close;
	}

	my $log = $self->lambda_bdg;
	$log =~ s/bdg$/out.txt/;
	my $command;

	# generate commands using bedtools and data2wig, which is faster than running
	# multiple instances of macs2 bdgcmp and bdgopt
	if ( $sfile and $lfile ) {

		# first step
		$command = sprintf
"%s unionbedg -header -names dlocal slocal llocal background -i %s %s %s %s > %s 2> $log ",
			$self->bedtools_app || 'bedtools',
			$dfile,
			$sfile,
			$lfile,
			$background_bdg,
			$self->sld_control_file;

		# second step
		$command .= sprintf
"&& %s --in %s --zero --fast --bdg --notrack --score 4-7 --method max --out %s ",
			$self->data2wig_app || 'data2wig.pl',
			$self->sld_control_file,
			$self->lambda_bdg;
		$command .= " 2>&1 >> $log ";

		# clean up
		$command .= sprintf "&& rm %s %s %s %s ",
			$sfile,
			$lfile,
			$background_bdg,
			$self->sld_control_file;
	}
	elsif ( $sfile and not $lfile ) {

		# first step
		$command = sprintf
"%s unionbedg -header -names dlocal slocal background -i %s %s %s > %s 2> $log ",
			$self->bedtools_app || 'bedtools',
			$dfile,
			$sfile,
			$background_bdg,
			$self->sld_control_file;

		# second step
		$command .= sprintf
"&& %s --in %s --zero --fast --bdg --notrack --score 4-6 --method max --out %s ",
			$self->data2wig_app || 'data2wig.pl',
			$self->sld_control_file,
			$self->lambda_bdg;
		$command .= " 2>&1 >> $log ";

		# clean up
		$command .= sprintf "&& rm %s %s %s ",
			$sfile,
			$background_bdg,
			$self->sld_control_file;
	}
	elsif ( not $sfile and $lfile ) {

		# first step
		$command = sprintf
"%s unionbedg -header -names dlocal llocal background -i %s %s %s > %s 2> $log ",
			$self->bedtools_app || 'bedtools',
			$dfile,
			$lfile,
			$background_bdg,
			$self->sld_control_file;

		# second step
		$command .= sprintf
"&& %s --in %s --zero --fast --bdg --notrack --score 4-6 --method max --out %s ",
			$self->data2wig_app || 'data2wig.pl',
			$self->sld_control_file,
			$self->lambda_bdg;
		$command .= " 2>&1 >> $log ";

		# clean up
		$command .= sprintf "&& rm %s %s %s ",
			$lfile,
			$background_bdg,
			$self->sld_control_file;
	}
	else {
		$self->crash(
			"programming error! how did we get here with no sfile and no lfile????");
	}

	$name2done->{ $self->lambda_bdg } = 1;
	return [ $command, $self->lambda_bdg, $log ];
}

sub convert_bw_to_bdg {
	my $self      = shift;
	my $name2done = shift;
	unless ( $self->bw2bdg_app =~ /\w+/ or $self->dryrun ) {
		croak "no bigWigToBedGraph application in path!\n";
	}
	my @commands;
	if ( $self->qvalue_bw and -e $self->qvalue_bw and not -e $self->qvalue_bdg ) {

		# the only way the q-value bigwig file is present and not the bedGraph is
		# if we're re-running a peak call and not a brand new pipeline
		# or somebody hasn't cleaned up somehow....
		my $log = $self->qvalue_bdg;
		$log =~ s/bdg$/bw2bdg.out.txt/;
		my $command = sprintf "%s %s %s 2>> $log",
			$self->bw2bdg_app || 'bigWigToBedGraph',
			$self->qvalue_bw,
			$self->qvalue_bdg;
		push @commands, [ $command, $self->qvalue_bdg, $log ];

		# there should not be any reason to continue in this situation
		return @commands;
	}
	if ( ( $self->chip_bw and -e $self->chip_bw and not -e $self->chip_bdg )
		or $self->dryrun )
	{
		# generate command if bigWig file exists or we're in dry run mode
		my $log = $self->chip_bdg;
		$log =~ s/bdg$/bw2bdg.out.txt/;
		my $command = sprintf
			"%s --in %s --out %s --skip \'%s\' --apply \'%s\' --bw2w %s 2>&1 > $log",
			$self->manwig_app || 'manipulate_wig.pl',
			$self->chip_bw,
			$self->chip_bdg,
			$self->chrskip,
			$self->chrskip,                          # use apply for tiny performance gain
			$self->bw2bdg_app || 'bigWigToBedGraph'; # to enforce bedGraph
		push @commands, [ $command, $self->chip_bdg, $log ];
	}
	if (
		(
				$self->lambda_bw
			and -e $self->lambda_bw
			and not -e $self->lambda_bdg
			and not exists $name2done->{ $self->lambda_bdg }
		)
		or $self->dryrun
		)
	{
		# generate command if bigWig exists and hasn't been done yet
		# or if we're in dry run mode
		my $log = $self->lambda_bdg;
		$log =~ s/bdg$/bw2bdg.out.txt/;
		my $command = sprintf "%s %s %s 2>> $log",
			$self->bw2bdg_app || 'bigWigToBedGraph',
			$self->lambda_bw,
			$self->lambda_bdg;
		push @commands, [ $command, $self->lambda_bdg, $log ];
		$name2done->{ $self->lambda_bdg } = 1;    # finished
	}
	return @commands;
}

sub generate_enrichment_commands {
	my $self = shift;
	unless ( $self->macs_app =~ /\w+/ or $self->dryrun ) {
		croak "no MACS2 application in path!\n";
	}
	my $chip_file   = $self->chip_bdg   || undef;
	my $lambda_file = $self->lambda_bdg || undef;
	unless ( $self->dryrun ) {
		if ( not $chip_file or not -e $chip_file ) {
			$self->crash("no ChIP fragment file $chip_file!\n");
		}
		if ( not $lambda_file or not -e $lambda_file ) {
			$self->crash("no control lambda fragment file $lambda_file!\n");
		}
	}

	# calculate scaling factor from RPM
	my $scale;
	if ( $self->samedepth ) {
		$scale = $self->targetdepth;
	}
	else {
		my @depths;
		if ( my $d = $self->seq_depth_for_file($chip_file) ) {
			push @depths, $d;
		}
		if ( my $d = $self->seq_depth_for_file($lambda_file) ) {
			push @depths, $d;
		}
		if (scalar @depths) {
			$scale = min(@depths);
		}
		else {
			$scale = $self->targetdepth;
		}
	}

	# generate command
	my $command = sprintf
		"%s bdgcmp -t %s -c %s -S %.4f -m qpois FE ",
		$self->macs_app || 'macs2',
		$chip_file,
		$lambda_file,
		$scale;
	if ( not $self->lambda and $lambda_file !~ /global_mean/ ) {
		$command .= "-p 1 ";    # add a pseudo count of 1 when doing explicit comparisons
	}
	$command .= sprintf "-o %s %s ", $self->qvalue_bdg, $self->fe_bdg;
	my $log = $self->qvalue_bdg;
	$log =~ s/qvalue \.bdg $/bdgcmp.out.txt/x;
	$command .= " 2> $log ";
	return [ $command, $self->qvalue_bdg, $log ];
}

sub generate_peakcall_commands {
	my $self = shift;
	unless ( $self->macs_app =~ /\w+/ or $self->dryrun ) {
		croak "no MACS2 application in path!\n";
	}

	my $qtrack = $self->qvalue_bdg;
	$self->crash("no qvalue bedGraph file defined!\n") if not $qtrack;
	unless ( $self->dryrun ) {
		$self->crash("no qvalue bedGraph file $qtrack!\n") if not -e $qtrack;
	}
	my $command = sprintf
		"%s bdgpeakcall -i %s -c %s -l %s -g %s --no-trackline -o %s ",
		$self->macs_app || 'macs2',
		$qtrack,
		$self->cutoff,
		$self->peaksize,
		$self->peakgap,
		$self->repmean_peak;
	my $log = $self->repmean_peak;
	$log =~ s/narrowPeak$/peakcall.out.txt/x;
	$command .= " 2> $log";

	if ( $self->broad ) {

		# sneak this in as an extra command
		# bdgbroadcall doesn't support writing without a trackline
		my $command2 = sprintf
			"%s bdgbroadcall -i %s -c %s -C %s -l %s -g %s -G %s -o %s ",
			$self->macs_app || 'macs2',
			$qtrack,
			$self->cutoff,
			$self->broadcut,
			$self->peaksize,
			$self->peakgap,
			$self->broadgap,
			$self->repmean_gappeak;
		my $log2 = $self->repmean_gappeak;
		$log2 =~ s/gappedPeak$/broadcall.out.txt/x;
		$command2 .= " 2> $log2";

		# return both narrow and broad commands
		return (
			[ $command,  $self->repmean_peak,    $log ],
			[ $command2, $self->repmean_gappeak, $log2 ]
		);
	}
	return [ $command, $self->repmean_peak, $log ];
}

sub generate_independent_peakcall_commands {
	my $self = shift;
	unless ( $self->macs_app =~ /\w+/ or $self->dryrun ) {
		croak "no MACS2 application in path!\n";
	}

	# ChIP bam files
	my $chip_number = scalar( $self->chip_use_bams );
	unless ( $chip_number >= 1 ) {
		$self->crash("no ChIP bam files available for independent peak call!");
	}

	# generic Macs2 options
	my $generic = q();
	if ( $self->paired ) {
		$generic .= '--format BAMPE ';
	}
	else {
		# single-end options
		# we're skipping modeling size - recommended to ensure consistent comparisons
		$generic .= sprintf "--format BAM --nomodel --extsize %s ", $self->fragsize;
		if ( $self->shiftsize ) {
			$generic .= sprintf "--shift %s ", $self->shiftsize;
		}
	}
	$generic .= '--keep-dup all ';    # we do our own de-deduplication
	if ( $self->lambda ) {
		$generic .= sprintf "--slocal %d --llocal %d ", $self->slocal, $self->llocal;
	}
	else {
		$generic .= '--nolambda ';
	}
	my $formatter = '%.' . sprintf( "%d", int( $self->cutoff + 0.5 ) ) . 'f';
	$generic .= sprintf
		"--qvalue $formatter --min-length %d --max-gap %d --gsize %d --outdir %s ",
		10**( -1 * $self->cutoff ),
		$self->peaksize,
		$self->peakgap,
		$self->genome,
		$self->dir;

	# broad specific options
	my $broad_generic = $generic;
	if ( $self->broad ) {
		$formatter = '%.' . sprintf( "%d", int( $self->broadcut + 1.5 ) ) . 'f';
		$broad_generic .= sprintf "--broad --broad-cutoff $formatter ",
			10**( -1 * $self->broadcut );
	}

	# generate commands for each ChIP replicate
	my @commands;
	if ( scalar( $self->control_use_bams ) == $chip_number ) {

		# perfect, equal numbers
		for ( my $i = 0; $i < $chip_number; $i++ ) {

			# narrow peak
			my $command = sprintf
				"%s callpeak --treatment %s --control %s --name %s $generic ",
				$self->macs_app || 'macs2',
				( $self->chip_use_bams )[$i],
				( $self->control_use_bams )[$i],
				( $self->chip_rep_names )[$i];
			my $out = ( $self->rep_peaks )[$i];
			my $log = $out;
			$log =~ s/narrowPeak$/narrowpeakcall.out.txt/x;
			$command .= "2> $log ";
			push @commands, [ $command, $out, $log ];

			# broad peaks
			if ( $self->broad ) {
				$command =
					sprintf
					"%s callpeak --treatment %s --control %s --name %s $broad_generic ",
					$self->macs_app || 'macs2',
					( $self->chip_use_bams )[$i],
					( $self->control_use_bams )[$i],
					( $self->chip_rep_names )[$i];
				$out = ( $self->rep_gappeaks )[$i];
				$log = $out;
				$log =~ s/gappedPeak$/broadcall.out.txt/x;
				$command .= "2> $log ";
				push @commands, [ $command, $out, $log ];
			}
		}
	}
	elsif ( scalar( $self->control_use_bams ) == 0 ) {

		# no control? ok
		for ( my $i = 0; $i < $chip_number; $i++ ) {

			# narrow peak
			my $command = sprintf
				"%s callpeak --treatment %s --name %s $generic ",
				$self->macs_app || 'macs2',
				( $self->chip_use_bams )[$i],
				( $self->chip_rep_names )[$i];
			my $out = ( $self->rep_peaks )[$i];
			my $log = $out;
			$log =~ s/narrowPeak$/narrowpeakcall.out.txt/x;
			$command .= "2> $log ";
			push @commands, [ $command, $out, $log ];

			# broad peak
			if ( $self->broad ) {
				$command = sprintf "%s callpeak --treatment %s --name %s $broad_generic ",
					$self->macs_app || 'macs2',
					( $self->chip_use_bams )[$i],
					( $self->chip_rep_names )[$i];
				$out = ( $self->rep_gappeaks )[$i];
				$log = $out;
				$log =~ s/gappedPeak$/broadcall.out.txt/x;
				$command .= "2> $log ";
				push @commands, [ $command, $out, $log ];
			}
		}
	}
	else {
		# use all the control bams for each separate chip
		for ( my $i = 0; $i < $chip_number; $i++ ) {

			# narrow peak
			my $command = sprintf
				"%s callpeak --treatment %s --control %s --name %s $generic ",
				$self->macs_app || 'macs2',
				( $self->chip_use_bams )[$i],
				join( q( ), ( $self->control_use_bams ) ),
				( $self->chip_rep_names )[$i];
			my $out = ( $self->rep_peaks )[$i];
			my $log = $out;
			$log =~ s/narrowPeak$/.narrowpeakcall.out.txt/x;
			$command .= "2> $log ";
			push @commands, [ $command, $out, $log ];

			# broad peak
			if ( $self->broad ) {
				$command =
					sprintf
					"%s callpeak --treatment %s --control %s --name %s $broad_generic ",
					$self->macs_app || 'macs2',
					( $self->chip_use_bams )[$i],
					join( q( ), ( $self->control_use_bams ) ),
					( $self->chip_rep_names )[$i];
				$out = ( $self->rep_gappeaks )[$i];
				$log = $out;
				$log =~ s/gappedPeak$/.broadcall.out.txt/x;
				$command .= "2> $log ";
				push @commands, [ $command, $out, $log ];
			}
		}
	}

	# finished
	return @commands;
}

sub generate_cleanpeak_commands {
	confess
"generate_cleanpeak_commands() is deprecated. See generate_peak_update_commands()";
}

sub generate_peak_update_commands {
	my $self = shift;
	unless ( $self->updatepeak_app =~ /\w+/ or $self->dryrun ) {
		croak "no update_peak_file.pl application in path!\n";
	}

	unless ( defined $self->repmean_peak ) {
		$self->crash("no peak file defined!");
	}

	# generate commands
	my @commands;
	my $update_options = sprintf "--enrich %s --qvalue %s --norm ",
		$self->logfe_bw, $self->qvalue_bw;

	# narrow peak
	my $log = $self->repmean_peak;
	$log =~ s/narrowPeak/update_peak.out.txt/;
	my $command = sprintf "%s --in %s %s --name %s --summit 2>&1 > $log",
		$self->updatepeak_app || 'update_peak_file.pl',
		$self->repmean_peak,
		$update_options,
		$self->job_name;
	my $summit = $self->repmean_peak;
	$summit =~ s/narrowPeak/summit.bed/;
	if ( $self->count_file_lines( $self->repmean_peak ) ) {
		push @commands, [ $command, $summit, $log ];
	}

	# broad peak
	if ($self->broad) {
		$log     = $self->repmean_gappeak . '.update_peak.out.txt';
		$command = sprintf "%s --in %s %s --name %s_broad 2>&1 > $log",
			$self->updatepeak_app || 'update_peak_file.pl',
			$self->repmean_gappeak,
			$update_options,
			$self->job_name;
		# since this command updates a pre-existing file without a new output file
		# we simply fake an output by indicating the stdout log file as the output
		if ( $self->count_file_lines( $self->repmean_gappeak ) ) {
			push @commands, [ $command, $summit, $log ];
		}
	}

	return @commands;
}

sub generate_independent_merge_peak_commands {
	my $self = shift;
	unless ( $self->intersect_app =~ /\w+/ or $self->dryrun ) {
		croak "no intersect_peaks.pl application in path!\n";
	}
	unless ( $self->bedtools_app =~ /\w+/ or $self->dryrun ) {
		croak "no bedtools application in path!\n";
	}

	# generate commands
	my @commands;
	if ( scalar( $self->rep_peaks ) == 1 ) {

		# only one replicate peak? nothing really to merge
		# to make it consistent with others, we will convert it to a simple bed file

		# narrow peak
		my $file  = ( $self->rep_peaks )[0];
		my $count = $self->count_file_lines($file);
		printf "  %s narrow peaks were called for %s replicate %s\n",
			format_with_commas($count), $self->job_name,
			( $self->chip_rep_names )[0];
		if ($count) {
			my $command = sprintf "%s --in %s --name %s --nosummit --norm --out %s ",
				$self->peak2bed_app || 'peak2bed.pl',
				$file,
				$self->job_name,
				$self->repmerge_peak;
			my $log = $file;
			$log =~ s/narrowPeak/peak2bed.out.txt/;
			$command .= " 2>&1 > $log";
			push @commands, [ $command, $self->repmerge_peak, $log ];
		}

		# broad peak
		if ( $self->broad ) {
			$file    = ( $self->rep_gappeaks )[0];
			$count = $self->count_file_lines($file);
			printf "  %s gapped peaks were called for %s replicate %s\n",
				format_with_commas($count), $self->job_name,
				( $self->chip_rep_names )[0];
			if ($count) {
				my $command = sprintf "%s --in %s --name %s --norm --out %s ",
					$self->peak2bed_app || 'peak2bed.pl',
					$file,
					$self->job_name,
					$self->repmerge_gappeak;
				my $log = $file;
				$log =~ s/gappedPeak/peak2bed.out.txt/;
				$command .= " 2>&1 > $log";
				push @commands, [ $command, $self->repmerge_gappeak, $log ];
			}
		}
	}
	elsif ( scalar( $self->rep_peaks ) > 1 ) {

		# more than one replicate peak, merge them
		# count the lines in each file, then merge non-empty files
		# narrow peaks first
		my @files;
		my @peaks = $self->rep_peaks;
		for my $i ( 0 .. $#peaks ) {
			my $count = $self->count_file_lines( $peaks[$i] );
			printf "  %s narrow peaks were called for %s replicate %s\n",
				format_with_commas($count), $self->job_name,
				( $self->chip_rep_names )[$i];
			if ($count) {
				push @files, $peaks[$i];
			}
		}
		if ( scalar @files == 1 ) {

			# only one replicate called peaks, use it
			# but only if there isn't a minimum number of overlap peaks
			if ( $self->minpeakover and $self->minpeakover > 1 ) {

				# can't do anything here
			}
			else {
				my $command = sprintf "%s --in %s --name %s --nosummit --norm --out %s ",
					$self->peak2bed_app || 'peak2bed.pl',
					$files[0],
					$self->job_name,
					$self->repmerge_peak;
				my $log = $files[0];
				$log =~ s/narrowPeak/peak2bed.out.txt/;
				$command .= " 2>&1 > $log";
				push @commands, [ $command, $self->repmerge_peak, $log ];
			}
		}
		elsif ( scalar @files > 1 ) {
			my $min = $self->minpeakover;
			unless ($min) {
				$min = scalar(@files) == 2 ? 2 : scalar(@files) - 1;
			}
			my $out = $self->repmerge_peak;
			$out =~ s/\.bed//;
			my $command = sprintf
				"%s --name %s --out %s --min %s --gap %s --bed %s --genome %s %s ",
				$self->intersect_app || 'intersect_peaks.pl',
				$self->job_name,
				$out,
				$min,
				$self->fragsize,
				$self->bedtools_app || 'bedtools',
				$self->chromofile,
				join( q( ), @files );
			my $log = sprintf "%s.intersect.out.txt", $out;
			$command .= " 2>&1 > $log ";
			push @commands, [ $command, $self->repmerge_peak, $log ];
		}

		# broad gapped peaks
		if ( $self->broad ) {
			@files = ();
			@peaks = $self->rep_gappeaks;
			foreach my $i ( 0 .. $#peaks ) {
				my $count = $self->count_file_lines( $peaks[$i] );
				printf "  %s gapped peaks were called for %s replicate %s\n",
					format_with_commas($count),
					$self->job_name,
					( $self->chip_rep_names )[$i];
				if ($count) {
					push @files, $peaks[$i];
				}
			}
			if ( scalar @files == 1 ) {

				# only one replicate called peaks, use it
				# but only if there isn't a minimum number of overlap peaks
				if ( $self->minpeakover and $self->minpeakover > 1 ) {

					# can't do anything here
				}
				else {
					my $command = sprintf "%s --in %s --name %s --norm --out %s ",
						$self->peak2bed_app || 'peak2bed.pl',
						$files[0],
						$self->job_name,
						$self->repmerge_gappeak;
					my $log = $files[0];
					$log =~ s/gappedPeak/peak2bed.out.txt/;
					$command .= " 2>&1 > $log";
					push @commands, [ $command, $self->repmerge_gappeak, $log ];
				}
			}
			elsif ( scalar @files > 1 ) {
				my $min = $self->minpeakover;
				unless ($min) {
					$min = scalar(@files) == 2 ? 2 : scalar(@files) - 1;
				}
				my $out = $self->repmerge_gappeak;
				$out =~ s/\.bed//;
				my $command = sprintf
					"%s --name %s --out %s --min %s --gap %s --bed %s --genome %s %s ",
					$self->intersect_app || 'intersect_peaks.pl',
					$self->job_name,
					$out,
					$min,
					$self->peaksize,
					$self->bedtools_app || 'bedtools',
					$self->chromofile,
					join( q( ), @files );
				my $log = $out . '.intersect.out.txt';
				$command .= " 2>&1 > $log ";
				push @commands, [ $command, $self->repmerge_gappeak, $log ];
			}
		}
	}
	return @commands;
}

sub generate_bdg2bw_commands {
	my $self      = shift;
	my $name2done = shift;
	unless ( $self->wig2bw_app =~ /\w+/ or $self->dryrun ) {
		croak "no wigToBigWig application in path!\n";
	}
	unless ( $self->manwig_app =~ /\w+/ or $self->dryrun ) {
		croak "no manipulate_wig.pl application in path!\n";
	}

	my @commands;

	# chip fragment files
	if ( $self->chip_bdg and $self->chip_bw ) {

		# we should have both files here
		if ( -e $self->chip_bdg and -e $self->chip_bw ) {

			# bigWig already exists so just delete the bdg
			push @commands, [ sprintf( "rm %s", $self->chip_bdg ), q(), q() ];
		}
		elsif ( ( -e $self->chip_bdg and not -e $self->chip_bw )
			or $self->dryrun )
		{
			# convert
			my $log = $self->chip_bw;
			$log =~ s/bw$/bdg2bw.out.txt/;
			my $command = sprintf
				"%s %s %s %s 2>&1 > $log && rm %s",
				$self->wig2bw_app || 'wigToBigWig',
				$self->chip_bdg,
				$self->chromofile,
				$self->chip_bw,
				$self->chip_bdg;
			push @commands, [ $command, $self->chip_bw, $log ];
		}
	}
	if (    $self->lambda_bdg
		and $self->lambda_bw
		and not exists $name2done->{ $self->lambda_bdg } )
	{
		if ( -e $self->lambda_bw ) {

			# we must have started with a lambda bigWig so remove the bedGraph
			push @commands, [ sprintf( "rm %s", $self->lambda_bdg ), q(), q() ];
			$name2done->{ $self->lambda_bdg } = 1;
		}
		else {
			# lambda control bigWig
			my $log = $self->lambda_bw;
			$log =~ s/bw$/bdg2bw.out.txt/;
			my $command = sprintf
				"%s %s %s %s 2>&1 > $log && rm %s",
				$self->wig2bw_app || 'wigToBigWig',
				$self->lambda_bdg,
				$self->chromofile,
				$self->lambda_bw,
				$self->lambda_bdg;
			push @commands, [ $command, $self->lambda_bw, $log ];

			# d control fragment bigWig
			if ( $self->d_control_bdg and $self->d_control_bw ) {
				$log = $self->d_control_bw;
				$log =~ s/bw$/bdg2bw.out.txt/;
				$command = sprintf
					"%s %s %s %s 2>&1 > $log && rm %s",
					$self->wig2bw_app || 'wigToBigWig',
					$self->d_control_bdg,
					$self->chromofile,
					$self->d_control_bw,
					$self->d_control_bdg;
				push @commands, [ $command, $self->d_control_bw, $log ];
			}
		}
		$name2done->{ $self->lambda_bdg } = 1;    # remember it's done
	}
	if ( $self->qvalue_bdg and $self->qvalue_bw ) {
		if ( -e $self->qvalue_bw ) {

			# qvalue bigWig already exists, probably because recalling peaks
			unless ( $self->savebdg ) {
				my $command = sprintf "rm %s", $self->qvalue_bdg;
				push @commands, [ $command, q(), q() ];
			}
		}
		elsif ( not -e $self->qvalue_bw or $self->dryrun ) {

			# convert bedGraph
			my $log = $self->qvalue_bw;
			$log =~ s/bw$/bdg2bw.out.txt/;
			my $command = sprintf
				"%s %s %s %s 2>&1 > $log ",
				$self->wig2bw_app || 'wigToBigWig',
				$self->qvalue_bdg,
				$self->chromofile,
				$self->qvalue_bw;

			unless ( $self->savebdg ) {
				$command .= sprintf " && rm %s", $self->qvalue_bdg;
			}
			push @commands, [ $command, $self->qvalue_bw, $log ];
		}
	}
	if ( $self->fe_bdg ) {

		# convert this to log2 Fold Enrichment because I like this better
		my $log = $self->logfe_bw;
		$log =~ s/bw$/bdg2bw.out.txt/;
		my $command = sprintf
			"%s --in %s --log 2 --place 4 --w2bw %s --chromo %s --out %s 2>&1 > $log ",
			$self->manwig_app || 'manipulate_wig.pl',
			$self->fe_bdg,
			$self->wig2bw_app || 'wigToBigWig',
			$self->chromofile,
			$self->logfe_bw;
		$command .= sprintf " && rm %s", $self->fe_bdg;
		push @commands, [ $command, $self->logfe_bw, $log ];
	}
	return @commands;
}

1;

=head1 NAME

Bio::MultiRepChIPSeq::Job - an object representing a ChIPSeq Experiment Job

=head1 DESCRIPTION

An object used to represent a ChIPSeq experiment consisting of one or more 
IP replicates and one or more corresponding reference control (Input) 
replicates as a Job in the Multi-replica ChIPSeq pipeline. Methods associated 
with this object represent in the processing steps for analyzing this experiment. 
Attributes include the input files and various intermediate and output files.

=head1 METHODS

=head2 Initialization

=over 4 

=item new()

Pass an array of seven items: 

=over 4

* A hash reference of options exported from a L<Bio::MultiRepChIPSeq::options> object 

* Name of the experiment

* Comma-delimited list of the ChIP bam files

* Comma-delimited list of the reference bam files, if any

* Comma-delimited list of the ChIP normalization scales, if any

* Comma-delimited list of the reference normalization scales, if any

* A chromosome-specific normalization scale, if defined

=back

=back

=head2 Array attributes

These methods return attributes consisting of an array, usually one or more filenames 
or values. Passing a value will add it to the list. 

=over 4

=item chip_bams

=item control_bams

=item chip_dedup_bams

=item chip_use_bams

=item control_dedup_bams

=item control_use_bams

=item chip_scale

=item control_scale

=back

=head2 Scalar attributes

These methods return scalar attributes of a single value. Passing a value will 
replace the existing one.

=over 4

=item job_name

=item chr_normfactor

=item chip_rep_names

=item rep_peaks

=item rep_gappeaks

=item chip_bw

=item chip_bdg

=item chip_count_bw

=item control_count_bw

=item d_control_bdg

=item d_control_bw

=item s_control_bdg

=item l_control_bdg

=item sld_control_file

=item lambda_bdg

=item lambda_bw

=item fe_bdg

=item logfe_bw

=item qvalue_bdg

=item qvalue_bw

=item repmean_peak

=item repmean_gappeak

=item repmerge_peak

=item repmerge_gappeak

=back

=head2 Action methods

These methods will generate commands to be acted on the files to process them. 
The commands can be executed by the L<Bio::MultiRepChIPSeq::Runner> object.

=over 4

=item generate_dedup_commands

=item find_dedup_bams

=item generate_bam2wig_frag_commands

=item generate_bam2wig_count_commands

=item generate_lambda_control_commands

=item convert_bw_to_bdg

=item generate_enrichment_commands

=item generate_peakcall_commands

=item generate_independent_peakcall_commands

=item generate_peak_update_commands

=item generate_independent_merge_peak_commands

=item generate_bdg2bw_commands

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

=head1 LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. 

