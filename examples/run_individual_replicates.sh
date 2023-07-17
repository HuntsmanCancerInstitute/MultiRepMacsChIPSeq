#!/bin/bash

# run one paired-end ChIP sample with multiple replicates individually

# cpu and job parameters are machine dependent and set low for this example.
# Due to extreme subsampling of alignments in example bam files, additional 
# parameters are supplied AND ARE NOT NORMALLY REQUIRED, including genome and plot options

# environment build paths – not needed in production
BLIB=${PWD}/../blib
export PATH=${BLIB}/script:$PATH
export PERL5LIB=${BLIB}/lib:$PERL5LIB

# clean up previous run results
rm -rf tup1_indiv

echo
echo "====================== Tup1 individual replicates ======================"
echo

multirep_macs2_pipeline.pl \
--chip data/Tup1_Ch1.bam \
--name Tup1_Ch1 \
--chip data/Tup1_Ch2.bam \
--name Tup1_Ch2 \
--chip data/Tup1_Ch3.bam \
--name Tup1_Ch3 \
--control data/Tup1_Input.bam \
--pe \
--dir tup1_indiv \
--out tup1_all \
--dupfrac 0.05 \
--optdist 2500 \
--size 275 \
--cutoff 3 \
--peaksize 200 \
--peakgap 100 \
--cpu 1 \
--job 4 \
--genome 230000 \
--plot \
--plot_frag 20000 \
--plot_qval 300


