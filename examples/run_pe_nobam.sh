#!/bin/bash

# example pipeline using only bigWig coverage files instead of bam files
# this is not preferred but allowed when coverage files are already available
# depending on how coverage tracks were generated this may or may not give
# desirable results

# the --targetdepth option must be set, either to 1 for unscaled coverage tracks
# or to the depth in millions for RPM scaled coverage tracks

# cpu and job parameters are machine dependent and set low for this example.
# Due to extreme subsampling of alignments in example bam files, additional 
# parameters are supplied AND ARE NOT NORMALLY REQUIRED, including genome options

# environment build paths – not needed in production
BLIB=${PWD}/../blib
export PATH=${BLIB}/script:$PATH
export PERL5LIB=${BLIB}/lib:$PERL5LIB

# clean up previous run results
rm -rf pe_nobam

# this requires that the run_pe.sh example has already been run
SRC=pe/Fragment
if [[ -e $SRC ]]
then
	echo
	echo "Found the bigWig coverage source files in $SRC"
else
	echo
	echo "Running the run_pe.sh pipeline first to generate source fragment files"
	./run_pe.sh
fi

echo
echo "====================== Peak call on bigWig only ======================"
echo

multirepchipseq.pl \
--chip $SRC/Rpd3.fragment.bw \
--control $SRC/Rpd3.lambda_control.bw \
--name Rpd3 \
--chip $SRC/Tup1.fragment.bw \
--control $SRC/Tup1.lambda_control.bw \
--name Tup1 \
--dir pe_nobam \
--out pe_all \
--targetdepth 0.1 \
--cutoff 3 \
--peaksize 200 \
--peakgap 100 \
--cpu 1 \
--job 4 \
--genome 230000 
