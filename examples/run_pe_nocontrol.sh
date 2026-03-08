#!/bin/bash

# two paired-end ChIP samples with no control
# as of v21 this uses large local lambda as default
# add --nolambda option for previous behavior

# cpu and job parameters are machine dependent and set low for this example.
# Due to extreme subsampling of alignments in example bam files, additional 
# parameters are supplied AND ARE NOT NORMALLY REQUIRED, including genome options

# environment build paths – not needed in production
BLIB=${PWD}/../blib
export PATH=${BLIB}/script:$PATH
export PERL5LIB=${BLIB}/lib:$PERL5LIB

# clean up previous run results
rm -rf pe_nocon

echo
echo "====================== paired-end with no control ======================"
echo

multirepchipseq.pl \
--chip data/Rpd3_Ch1.bam,data/Rpd3_Ch2.bam,data/Rpd3_Ch3.bam \
--name Rpd3 \
--chip data/Tup1_Ch1.bam,data/Tup1_Ch2.bam,data/Tup1_Ch3.bam \
--name Tup1 \
--pe \
--dir pe_nocon \
--out all_nocontrol \
--dupfrac 0.05 \
--optdist 2500 \
--cutoff 3 \
--peaksize 200 \
--peakgap 100 \
--cpu 1 \
--job 4 \
--genome 230000 

