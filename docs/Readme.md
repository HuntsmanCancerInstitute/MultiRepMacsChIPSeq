# MultiRepMacsChIPSeq

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Examples](Examples.md)|[Applications](applications.md)|[Install](Install.md)|

A multi-threaded wrapper for processing multi-replicate, multi-condition ChIPSeq samples

## Description

ChIPSeq and related methods, such as ATACSeq, is increasingly being used not just as
a means of discovery ("where is my factor binding?") but also as an assay for
experimental conditions ("how does my mutant affect factor X occupancy?"). Many
traditional ChIPSeq programs are not necessarily well equipped to handle multiple
replicates and/or sample conditions.

This is a wrapper application for processing ChIPSeq samples comprised of multiple
biological replicas and/or multiple conditions in a manner to make comparisons as
consistent and uniform as possible across samples.

## Rationale

The venerable [MACS2](https://pypi.org/project/MACS2/) application provides a robust
method of determining enrichment of ChIP fragments over input with a number of
advantages: fragment-based pileup of ChIP signal versus simple counts, single
base-pair resolution instead of sliding windows, estimation of local chromatin bias
using multiple window sizes, and more. However, Macs2 does not natively deal with
replicates, and comparing multiple conditions requires careful, manual execution 
of each set with identical parameters and/or complicated intersections. 

This package aims to automate Macs2 ChIPSeq peak calling with support for multiple
replicas and conditions while supporting newer normalization methods. Importantly, it
will output normalized, processed bigWig enrichment files for subsequent genic
analysis; numerous analytical, comparative, and quality contrrol metric plots; and
peak count tables ready for quantitative differential analysis. Finally, an overview
HTML report is generated with tables of numbers and selected QC and analytical plots.

**NOTE** This project currently uses MACS version 2. While
[MACS3](https://github.com/macs3-project/MACS) has been released, it has not yet been
evaluated and incorporated into this pipeline.

## Overview

See the [overview page](Overview.md) for a description of the steps involved in the 
pipeline.

## Installation

See the accompanying [Installation guide](INSTALL.md) for detailed notes on getting
the pipeline setup.

## Usage

See the [Usage Guide](Usage.md) for detailed notes on how to run the pipeline. For
variations in running the pipeline, such as with ATAC-Seq or Cut & Run, see the
[variations page](Variations.md).

A full [list of application menus](applications.md) is available. 

See the [examples page](Examples.md) for using the pipeline in common scenarios.

Users at HCI running the pipeline on local Linux servers can simply load the packages
into your environment using a `module` command. 

	module load multirepchipseq

This is also installed at CHPC at the University of Utah CHPC; please contact the
Cancer Bioinformatics shared resource for details.


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
