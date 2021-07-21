#!/usr/bin/perl

use strict;
use warnings;
use Module::Build;

my $build = Module::Build->new(
	dist_name         => 'MultiRepMacsChIPSeq',
	module_name       => 'MultiRepMacsChIPSeq',
	dist_version_from => 'lib/Bio/MultiRepChIPSeq/Runner.pm',
	dist_abstract     => 'Multiple-replica multiple-condition Macs2 ChIPSeq wrapper',
	dist_author       => 'Timothy Parnell <timothy.parnell@hci.utah.edu>',
	license           => 'artistic_2',
	meta_merge        => {
		resources     => {
			repository     => 'https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq'
		}
	},
	configure_requires => {
		'Module::Build'           => 0.42,
	},
	requires           => {
		'Bio::ToolBox'              => 1.69,
		'Bio::DB::HTS'              => 3.01,
		'Bio::DB::Big'              => 1,
		'Parallel::ForkManager'     => 1.19,
		'Set::IntervalTree'         => 0.11,
		'Set::IntSpan::Fast'        => 1.15,
		'Statistics::Descriptive'   => 3.0613,
	},
);

$build->create_build_script;
exit 0;
