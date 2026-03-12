package Bio::MultiRepChIPSeq;

use strict;
use Bio::MultiRepChIPSeq::Runner;

our $VERSION = 21.0;

sub new {
	return Bio::MultiRepChIPSeq::Runner->new();
}

1;

=head1 name

Bio::MultiRepChIPSeq - Multiple-replica multiple-condition ChIPSeq pipeline

=head1 SYNOPSIS

  use Bio::MultiRepChIPSeq;
  my $Runner = Bio::MultiRepChIPSeq->new();
  
  # set various parameters as appropriate
  # see Bio::MultiRepChIPSeq::options for a list of methods
  $Runner->paired(1);
  $Runner->mapq(30);
  
  # add first sample/condition job to the Runner
  $Runner->add_job(
     'MyChIP1',                            # name of experiment
     ['Rep1.bam', 'Rep2.bam', 'Rep3.bam'], # array of ChIP file names
     ['Input1.bam'],                       # array of reference file names
  );
  # add additional samples...
  
  # Run through pipeline steps
  $Runner->run_generate_chr_file();
  $Runner->run_dedup();
  $Runner->run_bam_check();
  $Runner->run_bam_fragment_conversion();
  $Runner->run_bam_count_conversion();
  $Runner->run_bw_conversion();
  $Runner->run_lambda_control();
  $Runner->run_bdgcmp();
  $Runner->run_call_peaks();
  $Runner->run_peak_merge();
  $Runner->run_bdg_conversion();
  $Runner->run_update_peaks();
  $Runner->run_rescore();
  $Runner->run_plot_peaks();


=head1 DESCRIPTION

This is a wrapper application for processing ChIPSeq samples comprised of
multiple biological replicas and/or multiple conditions in a manner to make
comparisons as consistent and uniform as possible across samples.

This package aims to automate Macs2 ChIPSeq peak calling with support for
multiple replicas and conditions while supporting newer normalization methods.
Importantly, it will output normalized, processed bigWig enrichment files for
subsequent genic analysis; numerous analytical, comparative, and quality
contrrol metric plots; and peak count tables ready for quantitative differential
analysis. Finally, an overview HTML report is generated with tables of numbers
and selected QC and analytical plots.

Typically, the L<multirepchipseq.pl> script, included in this distribution, is
used by end users to run the pipeline. The package includes numerous helper
scripts used during execution. Multiple external utilities are also required.

See L<https://huntsmancancerinstitute.github.io/MultiRepMacsChIPSeq> for details
on installation and usage.

=head1 METHODS

=head2 Initialization

=over 4

=item new

This returns a L<Bio::MultiRepChIPSeq::Runner> object, which is used to run the
pipeline.

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

=head1 LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. 

