#!/usr/bin/bash

################     Multi-replica ChIP-seq analysis with MACS2 ############

VERSION=2

################     USER VARIABLES     ####################################

# Modify the user variables

### Name
NAME=Experiment1


### Sorted indexed bam files
CHIPBAMFILES=( chip1.bam chip2.bam chip3.bam )
CONTROLBAMFILES=( input1.bam input2.bam input3.bam )

### Paired-end data
# set this to 1 for paired-end data
PAIRED=0


### Genome size
# human     2700000000
# mouse     1870000000
# zebrafish 1300000000
# fly       120000000
# c elegens 90000000
# yeast     12100000
GENOME=2700000000

### Fragment length
# should empirically check for each new dataset
FRAGMENTSIZE=300

### Target duplication rate
# set to 0 to skip duplicate removal
# set to 1 for duplicate subsample to target fraction
# set to 2 to remove all duplicates
REMOVE_DUPLICATES=1
MAXDUP=25
DUPFRACTION=0.1

### Peak calling parameters
QVAL=2
SLOCAL=1000
LLOCAL=10000
PEAKSIZE=300
PEAKGAP=200
TARGETDEPTH=25

### Chromosomes to skip
# this is a case sensitive regular expression and should be single double quoted
CHRSKIP='"chrM|MT|lambda|adapter|PhiX"'

### Black list bed file
# single doublequote this because it can be left null
BLACKLIST='""'


### Job number for parallel processing 
JOBNUMBER=3
THREADNUMBER=6


### Bin size in bp for making smaller bedgraphs
CHIPBIN=10
SLOCALBIN=50
LLOCALBIN=100



################    APPLICATION PATHS   ####################################

# BioToolBox programs are at https://github.com/tjparnell/biotoolbox
# extra scripts are at https://github.com/tjparnell/HCI-Scripts
# macs2 is at https://github.com/taoliu/MACS

### for local HCI servers
BAM2WIG=/home/BioApps/biotoolbox/bam2wig
MACS=/home/BioApps/macs2/macs2
MANWIG=/home/BioApps/biotoolbox/manipulate_wig.pl
BDG2BW=/tomato/dev/app/UCSC/bedGraphToBigWig
GNUPARALLEL=/usr/local/bin/parallel
PRINTCHR=/home/BioApps/biotoolbox/print_chromosome_lengths.pl
DEDUP=/home/BioApps/biotoolbox/bam_partial_dedup.pl




####################  MAIN SCRIPT   ########################################

echo
echo "========== Multi-replica Macs2 ChIPSeq Pipeline version $VERSION ============"
date
echo
echo "Name: $NAME"
echo "ChIP Bam files: ${CHIPBAMFILES[@]}"
echo "Control Bam files: ${CONTROLBAMFILES[@]}"
echo "Fragment size: $FRAGMENTSIZE bp"
if [ $REMOVE_DUPLICATES == 1 ]
then
	echo "Remove all duplicates"
elif [ $REMOVE_DUPLICATES == 2 ]
then
	echo "Subsample Duplicates"
	echo "Maximum duplicates allowed: $MAXDUP"
	echo "Subsample duplicates to fraction: $DUPFRACTION"
fi
echo "Skipped chromosomes: $CHRSKIP"
echo "Black list interval file: $BLACKLIST"
echo "Q-value threshold: $QVAL"
echo "Small lambda size: $SLOCAL bp"
echo "Large lambda size: $LLOCAL bp"
echo "Peak minimum size: $PEAKSIZE bp"
echo "Peak minimum gap size: $PEAKGAP bp"
echo
echo

# set parameter for paired end
PAIREDSTRING=''
if [ $PAIRED == 1 ]; then
	PAIREDSTRING=' --pe'
fi



#### Removing Duplicates
if [ $REMOVE_DUPLICATES == 1 ]
then
	echo "==== Subsampling duplicates to fraction $DUPFRACTION, maximum $MAXDUP"
	echo
	
	# run deduplication on bam files
	$GNUPARALLEL -j $JOBNUMBER -v -k \
	$DEDUP --max $MAXDUP --frac $DUPFRACTION --random \
	--seed 1 --in {} --out {.}.dedup.bam $PAIREDSTRING \
	':::' ${CHIPBAMFILES[@]} ${CONTROLBAMFILES[@]} 
	
	# check for new chip bam files
	NEWFILES=()
	for f in ${CHIPBAMFILES[@]}
	do
		if [ -s ${f%.bam}.dedup.bam ]
		then
			NEWFILES+=(${f%.bam}.dedup.bam)
		else
			# something went wrong, so use old one
			echo "PROBLEM: can't find ${f%.bam}.dedup.bam so using $f"
			NEWFILES+=($f)
		fi
	done
	CHIPBAMFILES=( ${NEWFILES[@]} )
	echo; echo "= NEW ChIP files ${CHIPBAMFILES[@]}"
	
	
	# check for new chip bam files
	NEWFILES=()
	for f in ${CONTROLBAMFILES[@]}
	do
		if [ -s ${f%.bam}.dedup.bam ]
		then
			NEWFILES+=(${f%.bam}.dedup.bam)
		else
			# something went wrong, so use old one
			echo "= PROBLEM: can't find ${f%.bam}.dedup.bam so using $f"
			NEWFILES+=($f)
		fi
	done
	CONTROLBAMFILES=( ${NEWFILES[@]} )
	echo; echo "= NEW Control files ${CONTROLBAMFILES[@]}"

elif [ $REMOVE_DUPLICATES == 2 ]
then
	echo "==== Removing all duplicates"
	echo
	
	# run deduplication on bam files
	$GNUPARALLEL -j $JOBNUMBER -v -k \
	$DEDUP --max 1 $PAIREDSTRING --in {} --out {.}.dedup.bam \
	':::' ${CHIPBAMFILES[@]} ${CONTROLBAMFILES[@]} 
	
	# check for new chip bam files
	NEWFILES=()
	for f in ${CHIPBAMFILES[@]}
	do
		if [ -e ${f%.bam}.dedup.bam ]
		then
			NEWFILES+=(${f%.bam}.dedup.bam)
		else
			# something went wrong, so use old one
			echo "PROBLEM: can't find ${f%.bam}.dedup.bam so using $f"
			NEWFILES+=($f)
		fi
	done
	CHIPBAMFILES=( ${NEWFILES[@]} )
	echo; echo "= NEW ChIP files ${CHIPBAMFILES[@]}"
	
	# check for new control bam files
	NEWFILES=()
	for f in ${CONTROLBAMFILES[@]}
	do
		if [ -e ${f%.bam}.dedup.bam ]
		then
			NEWFILES+=(${f%.bam}.dedup.bam)
		else
			# something went wrong, so use old one
			echo "PROBLEM: can't find ${f%.bam}.dedup.bam so using $f"
			NEWFILES+=($f)
		fi
	done
	CONTROLBAMFILES=( ${NEWFILES[@]} )
	echo; echo "= NEW Control files ${CONTROLBAMFILES[@]}"
fi



#### Generate Initial coverage tracks
# generation of lambda control based on 
# https://github.com/taoliu/MACS/wiki/Advanced%3A-Call-peaks-using-MACS2-subcommands

echo
echo "==== Generating Mean RPM ChIP coverage"

$BAM2WIG --cpu $THREADNUMBER --extend --extval $FRAGMENTSIZE \
--rpm --mean --bdg --bin $CHIPBIN --chrskip $CHRSKIP --blacklist $BLACKLIST \
$PAIREDSTRING --out $NAME.extend.bdg ${CHIPBAMFILES[@]}

echo
echo "==== Generating Mean RPM d control"

$BAM2WIG --cpu $THREADNUMBER --cspan --extval $FRAGMENTSIZE \
--rpm --mean --bdg --bin $CHIPBIN --chrskip $CHRSKIP --blacklist $BLACKLIST \
$PAIREDSTRING --out $NAME.dlocal.bdg ${CONTROLBAMFILES[@]}

echo
echo "==== Generating Mean RPM slocal control"

SSCALE=$(echo "$FRAGMENTSIZE / $SLOCAL" | bc -l)
echo " Small local lambda scale $SSCALE"
$BAM2WIG --cpu $THREADNUMBER --cspan --extval $SLOCAL --rpm \
--mean --bdg --bin $SLOCALBIN --chrskip $CHRSKIP --blacklist $BLACKLIST \
--scale $SSCALE $PAIREDSTRING --out $NAME.slocal.bdg ${CONTROLBAMFILES[@]}

echo
echo "==== Generating Mean RPM llocal control"

LSCALE=$(echo "$FRAGMENTSIZE / $LLOCAL" | bc -l)
echo " Large local lambda scale $LSCALE"
$BAM2WIG --cpu $THREADNUMBER --cspan --extval $LLOCAL --rpm \
--mean --bdg --bin $LLOCALBIN --chrskip $CHRSKIP --blacklist $BLACKLIST \
--scale $LSCALE $PAIREDSTRING --out $NAME.llocal.bdg ${CONTROLBAMFILES[@]}

echo
echo "==== Generating lambda control"

$MACS bdgcmp -m max -t $NAME.slocal.bdg -c $NAME.llocal.bdg -o $NAME.sllocal.bdg \
&& rm $NAME.slocal.bdg $NAME.llocal.bdg

$MACS bdgcmp -m max -t $NAME.sllocal.bdg -c $NAME.dlocal.bdg -o $NAME.sldlocal.bdg \
&& rm $NAME.sllocal.bdg $NAME.dlocal.bdg

# genomic background based on 1 million reads
BACKGROUND=$(expr 1000000 \* $FRAGMENTSIZE)
BACKGROUND=$(echo "$BACKGROUND / $GENOME" | bc -l)
echo "= Genome background depth is $BACKGROUND"

$MACS bdgopt -i $NAME.sldlocal.bdg -m max -p $BACKGROUND -o $NAME.control_lambda.bdg \
&& rm $NAME.sldlocal.bdg





#### Generate Initial coverage tracks
echo
echo "==== Generating enrichment tracks"

$MACS bdgcmp -t $NAME.extend.bdg -c $NAME.control_lambda.bdg -S $TARGETDEPTH \
-m qpois FE -o $NAME.qvalue.bdg $NAME.FE.bdg


echo
echo "==== Converting log2 fold enrichment tracks"

$MANWIG --in $NAME.FE.bdg --log 2 --place 4 --out $NAME.log2FE.bdg \
&& rm $NAME.FE.bdg




#### Call peaks
echo
echo "==== Calling peaks"
$MACS bdgpeakcall -i $NAME.qvalue.bdg -c $QVAL -l $PEAKSIZE -g $PEAKGAP \
--no-trackline -o $NAME.narrowPeak



### Convert 
echo
echo "==== Converting tracks to bigWig"
$PRINTCHR ${CHIPBAMFILES[0]} > chrom.sizes

$GNUPARALLEL -k -v -j $JOBNUMBER \
$BDG2BW {} chrom.sizes {.}.bw \
':::' $NAME.extend.bdg $NAME.control_lambda.bdg $NAME.qvalue.bdg $NAME.log2FE.bdg

rm chrom.sizes

echo
date
echo "====== finished ======="

# Timothy J. Parnell, PhD
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112







