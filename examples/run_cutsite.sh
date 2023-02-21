#!/bin/bash

# two paired-end samples run with ATAC-Seq cut site analysis

# environment build paths – not needed in production
BLIB=${PWD}/../blib
export PATH=${BLIB}/script:$PATH
export PERL5LIB=${BLIB}/lib:$PERL5LIB

# clean up previous run results
rm -rf cutsite

echo
echo "====================== paired-end cut site ======================"
echo

multirep_macs2_pipeline.pl \
--chip data/Rpd3_Ch1.bam,data/Rpd3_Ch2.bam,data/Rpd3_Ch3.bam \
--name Rpd3 \
--chip data/Tup1_Ch1.bam,data/Tup1_Ch2.bam,data/Tup1_Ch3.bam \
--name Tup1 \
--dir cutsite \
--out all_cut \
--genome 230000 \
--deduppair \
--dupfrac 0.1 \
--optdist 2500 \
--shift -50 \
--size 100 \
--cutoff 3 \
--peaksize 150 \
--peakgap 50 \
--plot \
--cpu 1 \
--job 4 
