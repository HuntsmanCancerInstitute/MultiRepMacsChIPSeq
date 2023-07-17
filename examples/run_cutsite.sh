#!/bin/bash

# two paired-end samples run with ATAC-Seq cut site analysis
# use new atac option introduced in v18.1 to generalize parameters

# cpu and job parameters are machine dependent and set low for this example.
# Due to extreme subsampling of alignments in example bam files, additional 
# parameters are supplied AND ARE NOT NORMALLY REQUIRED, including genome and plot options

# environment build paths – not needed in production
BLIB=${PWD}/../blib
export PATH=${BLIB}/script:$PATH
export PERL5LIB=${BLIB}/lib:$PERL5LIB

# clean up previous run results
rm -rf cutsite

echo
echo "====================== ATAC-Seq cut site ======================"
echo

multirep_macs2_pipeline.pl \
--chip data/Rpd3_Ch1.bam,data/Rpd3_Ch2.bam,data/Rpd3_Ch3.bam \
--name Rpd3 \
--chip data/Tup1_Ch1.bam,data/Tup1_Ch2.bam,data/Tup1_Ch3.bam \
--name Tup1 \
--dir cutsite \
--out all_cut \
--dupfrac 0.1 \
--optdist 2500 \
--atac \
--cutoff 3 \
--cpu 1 \
--job 4 \
--genome 230000 \
--plot \
--plot_frag 25000 \
--plot_qval 300

