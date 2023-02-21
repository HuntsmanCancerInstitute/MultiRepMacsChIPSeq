#!/bin/bash

# run two paired-end ChIP samples, with independent peak calls for the one with replicates

# environment build paths – not needed in production
BLIB=${PWD}/../blib
export PATH=${BLIB}/script:$PATH
export PERL5LIB=${BLIB}/lib:$PERL5LIB

# clean up previous run results
rm -rf indep_pe

echo
echo "====================== independent paired-end ======================"
echo

multirep_macs2_pipeline.pl \
--chip data/Rpd3_Ch1.bam,data/Rpd3_Ch2.bam,data/Rpd3_Ch3.bam \
--control data/Rpd3_Input.bam \
--name Rpd3 \
--chip data/Tup1_Ch1.bam \
--control data/Tup1_Input.bam \
--name Tup1 \
--pe \
--dir indep_pe \
--out pe_all \
--dupfrac 0.05 \
--optdist 2500 \
--size 275 \
--cutoff 3 \
--independent \
--peaksize 200 \
--peakgap 100 \
--plot \
--cpu 1 \
--job 4 

