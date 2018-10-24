#!/usr/bin/perl

use strict;
use warnings;
use Module::Build;

my $build = Module::Build->new(
	dist_name         => 'MultiRepMacsChIPSeq',
	module_name       => 'MultiRepMacsChIPSeq',
	dist_version_from => 'bin/multirep_macs2_pipeline.pl',
	dist_abstract     => 'Multiple-replica multiple-condition Macs2 ChIPSeq wrapper',
	dist_author       => 'Timothy Parnell <timothy.parnell@hci.utah.edu>',
	license           => 'artistic_2',
	meta_merge        => {
		resources     => {
			repository     => 'https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq'
		}
	},
	configure_requires => {
		'Module::Build'           => 0,
	},
	requires           => {
		'Bio::ToolBox'              => 1.63,
		'Parallel::ForkManager'     => 0,
		'Bio::DB::HTS'              => 0,
		'Bio::DB::Big'              => 0,
	},
);

$build->create_build_script;
exit 0;
