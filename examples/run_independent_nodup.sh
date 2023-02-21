#!/bin/bash

# run two paired-end ChIP samples, with independent peak calls, without duplication

# environment build paths – not needed in production
BLIB=${PWD}/../blib
export PATH=${BLIB}/script:$PATH
export PERL5LIB=${BLIB}/lib:$PERL5LIB

# clean up previous run results
rm -rf nodup_pe

echo
echo "================== Independent Paired-end without duplication =================="
echo

multirep_macs2_pipeline.pl \
--chip data/Rpd3_Ch1.bam,data/Rpd3_Ch2.bam,data/Rpd3_Ch3.bam \
--control data/Rpd3_Input.bam \
--name Rpd3 \
--chip data/Tup1_Ch1.bam,data/Tup1_Ch2.bam,data/Tup1_Ch3.bam \
--control data/Tup1_Input.bam \
--name Tup1 \
--pe \
--dir nodup_pe \
--out all_nodup \
--genome 230000 \
--nodedup \
--independent \
--cutoff 3 \
--peaksize 200 \
--peakgap 100 \
--plot \
--cpu 1 \
--job 4 


