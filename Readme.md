# Multi-Replica Macs ChIPSeq Wrapper

A multi-threaded wrapper for processing multi-replicate, multi-condition ChIPSeq samples

## Description

ChIPSeq and related methods, such as ATACSeq and Cut-and-Run, is increasingly being
used not just as a means of discovery ("where is my factor binding?") but also as an
assay for experimental conditions ("how does my mutant affect factor X occupancy?").
Many traditional ChIPSeq programs are not necessarily well equipped to handle
multiple replicates and/or sample conditions.

This is a wrapper application for processing ChIPSeq samples comprised of multiple
biological replicas and/or multiple conditions in a manner to make comparisons as
consistent and uniform as possible across samples. 

## Documentation

For full documentation, including overview, usage guide, examples, and installation
guide, see the [Documentation pages](https://huntsmancancerinstitute.github.io/MultiRepMacsChIPSeq).

## Installation

This is a Perl package and may be built using the standard incantation.

	perl Build.PL
	./Build
	./Build install

For detailed notes on accessory packages and dependencies, see the 
[Installation guide](https://huntsmancancerinstitute.github.io/MultiRepMacsChIPSeq/Install.hmtl).

## AUTHOR

    Timothy J. Parnell, PhD
    Cancer Bioinformatics Shared Resource
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

## LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. For details, see the
full text of the license in the file LICENSE.

This package is distributed in the hope that it will be useful, but it
is provided "as is" and without any express or implied warranties. For
details, see the full text of the license in the file LICENSE.
