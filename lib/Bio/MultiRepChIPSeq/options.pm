package Bio::MultiRepChIPSeq::options;

use strict;
use Carp;
use File::Which;

our $VERSION = 18.0;

sub init_options {
	my $class = shift;
	my %opts  = (
		dir         => './',
		in          => undef,
		out         => 'merged',
		name        => [],
		chip        => [],
		control     => [],
		chscale     => [],
		coscale     => [],
		genome      => 0,
		species     => q(),
		chromofile  => undef,
		mapq        => 0,
		paired      => 0,
		fraction    => 0,
		minsize     => 50,
		maxsize     => 500,
		chrskip     => "chrM|MT|alt|Adapter|Lambda|PhiX",
		blacklist   => undef,
		dedup       => 1,
		maxdup      => undef,                               # old option
		maxdepth    => undef,
		dupfrac     => 0.05,
		optdist     => 0,
		deduppair   => 0,
		fragsize    => 250,
		shiftsize   => 0,
		slocal      => 1000,
		llocal      => 10000,
		cutoff      => 2,
		targetdep   => undef,
		peaksize    => undef,
		peakgap     => undef,
		broad       => 0,
		broadcut    => 0.5,
		broadgap    => undef,
		gaplink     => undef,
		lambda      => 1,
		independent => 0,
		cpu         => 4,
		job         => 2,
		chipbin     => 10,
		slocalbin   => 50,
		llocalbin   => 100,
		chrapply    => undef,
		chrnorm     => [],
		rawcounts   => 0,
		savebam     => 0,
		savebdg     => 0,
		binsize     => 40,
		genomewin   => 0,
		discard     => 10,
		repmean     => 0,
		plot        => 1,
		dryrun      => 0,
		organize    => 1,
		bam2wig     => sprintf( "%s", which 'bam2wig.pl' ),
		bamdedup    => sprintf( "%s", which 'bam_partial_dedup.pl' ),
		macs        => sprintf( "%s", which 'macs2' ),
		manwig      => sprintf( "%s", which 'manipulate_wig.pl' ),
		wig2bw      => sprintf( "%s", which 'wigToBigWig' ),
		bw2bdg      => sprintf( "%s", which 'bigWigToBedGraph' ),
		bedtools    => sprintf( "%s", which 'bedtools' ),
		getdata     => sprintf( "%s", which 'get_datasets.pl' ),
		getrel      => sprintf( "%s", which 'get_relative_data.pl' ),
		geteff      => sprintf( "%s", which 'get_chip_efficiency.pl' ),
		printchr    => sprintf( "%s", which 'print_chromosome_lengths.pl' ),
		data2wig    => sprintf( "%s", which 'data2wig.pl' ),
		meanbdg     => sprintf( "%s", which 'generate_mean_bedGraph.pl' ),
		intersect   => sprintf( "%s", which 'intersect_peaks.pl' ),
		peak2bed    => sprintf( "%s", which 'peak2bed.pl' ),
		combrep     => sprintf( "%s", which 'combine_replicate_data.pl' ),
		plotpeak    => sprintf( "%s", which 'plot_peak_figures.R' ),
		rscript     => sprintf( "%s", which 'Rscript' ),
		reportmap   => sprintf( "%s", which 'report_mappable_space.pl' ),
		samtools    => sprintf( "%s", which 'samtools' ),
		help        => 0,
	);
	return \%opts;
}

sub options {
	my $self = shift;
	return $self->{opts};
}

# sub xxx {
# 	my $self = shift;
# 	$self->{opts}{xxx} = shift if @_;
# 	return $self->{opts}{xxx};
# }

sub dir {
	my $self = shift;
	$self->{opts}{dir} = shift if @_;
	return $self->{opts}{dir};
}

sub in {
	my $self = shift;
	$self->{opts}{in} = shift if @_;
	return $self->{opts}{in};
}

sub out {
	my $self = shift;
	$self->{opts}{out} = shift if @_;
	return $self->{opts}{out};
}

sub name {
	my $self = shift;
	push @{ $self->{opts}{name} }, @_ if @_;
	return @{ $self->{opts}{name} };
}

sub chip {
	my $self = shift;
	push @{ $self->{opts}{chip} }, @_ if @_;
	return @{ $self->{opts}{chip} };
}

sub control {
	my $self = shift;
	push @{ $self->{opts}{control} }, @_ if @_;
	return @{ $self->{opts}{control} };
}

sub chscale {
	my $self = shift;
	push @{ $self->{opts}{chscale} }, @_ if @_;
	return @{ $self->{opts}{chscale} };
}

sub coscale {
	my $self = shift;
	push @{ $self->{opts}{coscale} }, @_ if @_;
	return @{ $self->{opts}{coscale} };
}

sub genome {
	my $self = shift;
	$self->{opts}{genome} = shift if @_;
	return $self->{opts}{genome};
}

sub species {
	my $self = shift;
	$self->{opts}{species} = shift if @_;
	return $self->{opts}{species};
}

sub chromofile {
	my $self = shift;
	$self->{opts}{chromofile} = shift if @_;
	return $self->{opts}{chromofile};
}

sub paired {
	my $self = shift;
	$self->{opts}{paired} = shift if @_;
	return $self->{opts}{paired};
}

sub mapq {
	my $self = shift;
	$self->{opts}{mapq} = shift if @_;
	return $self->{opts}{mapq};
}

sub fraction {
	my $self = shift;
	$self->{opts}{fraction} = shift if @_;
	return $self->{opts}{fraction};
}

sub minsize {
	my $self = shift;
	$self->{opts}{minsize} = shift if @_;
	return $self->{opts}{minsize};
}

sub maxsize {
	my $self = shift;
	$self->{opts}{maxsize} = shift if @_;
	return $self->{opts}{maxsize};
}

sub dedup {
	my $self = shift;
	$self->{opts}{dedup} = shift if @_;
	return $self->{opts}{dedup};
}

sub maxdup {
	my $self = shift;
	$self->{opts}{maxdup} = shift if @_;
	return $self->{opts}{maxdup};
}

sub maxdepth {
	my $self = shift;
	$self->{opts}{maxdepth} = shift if @_;
	return $self->{opts}{maxdepth};
}

sub dupfrac {
	my $self = shift;
	$self->{opts}{dupfrac} = shift if @_;
	return $self->{opts}{dupfrac};
}

sub optdist {
	my $self = shift;
	$self->{opts}{optdist} = shift if @_;
	return $self->{opts}{optdist};
}

sub deduppair {
	my $self = shift;
	$self->{opts}{deduppair} = shift if @_;
	return $self->{opts}{deduppair};
}

sub savebam {
	my $self = shift;
	$self->{opts}{savebam} = shift if @_;
	return $self->{opts}{savebam};
}

sub fragsize {
	my $self = shift;
	$self->{opts}{fragsize} = shift if @_;
	return $self->{opts}{fragsize};
}

sub shiftsize {
	my $self = shift;
	$self->{opts}{shiftsize} = shift if @_;
	return $self->{opts}{shiftsize};
}

sub slocal {
	my $self = shift;
	$self->{opts}{slocal} = shift if @_;
	return $self->{opts}{slocal};
}

sub llocal {
	my $self = shift;
	$self->{opts}{llocal} = shift if @_;
	return $self->{opts}{llocal};
}

sub cutoff {
	my $self = shift;
	$self->{opts}{cutoff} = shift if @_;
	return $self->{opts}{cutoff};
}

sub targetdep {
	my $self = shift;
	$self->{opts}{targetdep} = shift if @_;
	return $self->{opts}{targetdep};
}

sub peaksize {
	my $self = shift;
	$self->{opts}{peaksize} = shift if @_;
	return $self->{opts}{peaksize};
}

sub peakgap {
	my $self = shift;
	$self->{opts}{peakgap} = shift if @_;
	return $self->{opts}{peakgap};
}

sub broad {
	my $self = shift;
	$self->{opts}{broad} = shift if @_;
	return $self->{opts}{broad};
}

sub broadcut {
	my $self = shift;
	$self->{opts}{broadcut} = shift if @_;
	return $self->{opts}{broadcut};
}

sub broadgap {
	my $self = shift;
	$self->{opts}{broadgap} = shift if @_;
	return $self->{opts}{broadgap};
}

sub gaplink {
	my $self = shift;
	$self->{opts}{gaplink} = shift if @_;
	return $self->{opts}{gaplink};
}

sub lambda {
	my $self = shift;
	$self->{opts}{lambda} = shift if @_;
	return $self->{opts}{lambda};
}

sub independent {
	my $self = shift;
	$self->{opts}{independent} = shift if @_;
	return $self->{opts}{independent};
}

sub chrskip {
	my $self = shift;
	$self->{opts}{chrskip} = shift if @_;
	return $self->{opts}{chrskip};
}

sub blacklist {
	my $self = shift;
	$self->{opts}{blacklist} = shift if @_;
	return $self->{opts}{blacklist};
}

sub cpu {
	my $self = shift;
	$self->{opts}{cpu} = shift if @_;
	return $self->{opts}{cpu};
}

sub job {
	my $self = shift;
	$self->{opts}{job} = shift if @_;
	return $self->{opts}{job};
}

sub chipbin {
	my $self = shift;
	$self->{opts}{chipbin} = shift if @_;
	return $self->{opts}{chipbin};
}

sub slocalbin {
	my $self = shift;
	$self->{opts}{slocalbin} = shift if @_;
	return $self->{opts}{slocalbin};
}

sub llocalbin {
	my $self = shift;
	$self->{opts}{llocalbin} = shift if @_;
	return $self->{opts}{llocalbin};
}

sub chrnorm {
	my $self = shift;
	push @{ $self->{opts}{chrnorm} }, @_ if @_;
	return @{ $self->{opts}{chrnorm} };
}

sub chrapply {
	my $self = shift;
	$self->{opts}{chrapply} = shift if @_;
	return $self->{opts}{chrapply};
}

sub rawcounts {
	my $self = shift;
	$self->{opts}{rawcounts} = shift if @_;
	return $self->{opts}{rawcounts};
}

sub savebdg {
	my $self = shift;
	$self->{opts}{savebdg} = shift if @_;
	return $self->{opts}{savebdg};
}

sub binsize {
	my $self = shift;
	$self->{opts}{binsize} = shift if @_;
	return $self->{opts}{binsize};
}

sub genomewin {
	my $self = shift;
	$self->{opts}{genomewin} = shift if @_;
	return $self->{opts}{genomewin};
}

sub discard {
	my $self = shift;
	$self->{opts}{discard} = shift if @_;
	return $self->{opts}{discard};
}

sub repmean {
	my $self = shift;
	$self->{opts}{repmean} = shift if @_;
	return $self->{opts}{repmean};
}

sub plot {
	my $self = shift;
	$self->{opts}{plot} = shift if @_;
	return $self->{opts}{plot};
}

sub dryrun {
	my $self = shift;
	$self->{opts}{dryrun} = shift if @_;
	return $self->{opts}{dryrun};
}

sub organize {
	my $self = shift;
	$self->{opts}{organize} = shift if @_;
	return $self->{opts}{organize};
}

sub bam2wig_app {
	my $self = shift;
	$self->{opts}{bam2wig} = shift if @_;
	return $self->{opts}{bam2wig};
}

sub bamdedup_app {
	my $self = shift;
	$self->{opts}{bamdedup} = shift if @_;
	return $self->{opts}{bamdedup};
}

sub macs_app {
	my $self = shift;
	$self->{opts}{macs} = shift if @_;
	return $self->{opts}{macs};
}

sub manwig_app {
	my $self = shift;
	$self->{opts}{manwig} = shift if @_;
	return $self->{opts}{manwig};
}

sub wig2bw_app {
	my $self = shift;
	$self->{opts}{wig2bw} = shift if @_;
	return $self->{opts}{wig2bw};
}

sub bw2bdg_app {
	my $self = shift;
	$self->{opts}{bw2bdg} = shift if @_;
	return $self->{opts}{bw2bdg};
}

sub bedtools_app {
	my $self = shift;
	$self->{opts}{bedtools} = shift if @_;
	return $self->{opts}{bedtools};
}

sub getdata_app {
	my $self = shift;
	$self->{opts}{getdata} = shift if @_;
	return $self->{opts}{getdata};
}

sub getrel_app {
	my $self = shift;
	$self->{opts}{getrel} = shift if @_;
	return $self->{opts}{getrel};
}

sub geteff_app {
	my $self = shift;
	$self->{opts}{geteff} = shift if @_;
	return $self->{opts}{geteff};
}

sub printchr_app {
	my $self = shift;
	$self->{opts}{printchr} = shift if @_;
	return $self->{opts}{printchr};
}

sub data2wig_app {
	my $self = shift;
	$self->{opts}{data2wig} = shift if @_;
	return $self->{opts}{data2wig};
}

sub meanbdg_app {
	my $self = shift;
	$self->{opts}{meanbdg} = shift if @_;
	return $self->{opts}{meanbdg};
}

sub intersect_app {
	my $self = shift;
	$self->{opts}{intersect} = shift if @_;
	return $self->{opts}{intersect};
}

sub peak2bed_app {
	my $self = shift;
	$self->{opts}{peak2bed} = shift if @_;
	return $self->{opts}{peak2bed};
}

sub combrep_app {
	my $self = shift;
	$self->{opts}{combrep} = shift if @_;
	return $self->{opts}{combrep};
}

sub plotpeak_app {
	my $self = shift;
	$self->{opts}{plotpeak} = shift if @_;
	return $self->{opts}{plotpeak};
}

sub rscript_app {
	my $self = shift;
	$self->{opts}{rscript} = shift if @_;
	return $self->{opts}{rscript};
}

sub reportmap_app {
	my $self = shift;
	$self->{opts}{reportmap} = shift if @_;
	return $self->{opts}{reportmap};
}

sub samtools_app {
	my $self = shift;
	$self->{opts}{samtools} = shift if @_;
	return $self->{opts}{samtools};
}

1;

=head1 NAME

Bio::MultiRepChIPSeq::options - an object to hold the pipeline options

=head1 DESCRIPTION

This is a base module inherited by the others for holding general configuration 
parameters and options. The values can be exported as a hash reference for use in
importing command line options with L<Getopt::Long>. The values can be get or set 
using methods generally named after the hash key (with the exception of application 
paths). 

This is not used directly. 

=head1 METHODS

=head2 Initialization

=over 4

=item init_options

Initialize the options hash has and identify application paths in the current 
environment C<PATH>. 

=item options

Export the hash as a reference.

=back

=head2 Attributes

=over 4

=item dir

=item in

=item out

=item name

=item chip

=item control

=item chscale

=item coscale

=item genome

=item species

=item chromofile

=item paired

=item mapq

=item fraction

=item minsize

=item maxsize

=item dedup

=item maxdup

=item maxdepth

=item dupfrac

=item optdist

=item deduppair

=item savebam

=item fragsize

=item shiftsize

=item slocal

=item llocal

=item cutoff

=item targetdep

=item peaksize

=item peakgap

=item broad

=item broadcut

=item broadgap

=item gaplink

=item lambda

=item independent

=item chrskip

=item blacklist

=item cpu

=item job

=item chipbin

=item slocalbin

=item llocalbin

=item chrnorm

=item chrapply

=item rawcounts

=item savebdg

=item binsize

=item genomewin

=item discard

=item repmean

=item plot

=item dryrun

=item organize

=item bam2wig_app

=item bamdedup_app

=item macs_app

=item manwig_app

=item wig2bw_app

=item bw2bdg_app

=item bedtools_app

=item getdata_app

=item getrel_app

=item geteff_app

=item printchr_app

=item data2wig_app

=item meanbdg_app

=item intersect_app

=item peak2bed_app

=item combrep_app

=item plotpeak_app

=item rscript_app

=item reportmap_app

=item samtools_app

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

=head1 LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. 


