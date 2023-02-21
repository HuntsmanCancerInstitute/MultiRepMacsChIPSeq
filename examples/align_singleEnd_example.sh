#!/bin/bash

# for HCI pysano purposes only - delete this and preceding lines and rename as cmd.txt
#e your.name@hci.utah.edu
#c kingspeak
#t 12

# single-end ChIPSeq
# alignment with Novoalign followed by duplication metrics and fragment size estimation

# set the base name of the output
NAME=MYNAME

# mouse Novoalign index
DATA=/tomato/dev/data
INDEX=$DATA/Mouse/Mm10/mm10.standard.nov.illumina.nix
SIZE=2500000000

# optical distance
# use 100 for HiSeq runs or 2500 for NovaSeq runs
OPTDISTANCE=100

# application paths
NOVO_APP=/tomato/dev/app/novoalign/4.03.01/novoalign
SAM_APP=/tomato/dev/app/samtools/1.16/samtools
DEDUP_APP=/tomato/dev/app/modulesoftware/bam_partial_dedup
BAMWIG_APP=/tomato/dev/app/modulesoftware/bam2wig
MACS_APP=/tomato/dev/app/modulesoftware/macs2

# align
echo "=== Aligning ==="
$NOVO_APP \
-d $INDEX \
-o SAM \
-a $ADAPTF $ADAPTR \
--tune HiSeq \
-f *.gz | \
$SAM_APP view -b -o $NAME.raw.bam -

# sort and index
echo "=== Sorting ==="
$SAM_APP sort -m 4G -@ 8 -o $NAME.bam $NAME.raw.bam \
&& $SAM_APP index -@ 8 $NAME.bam \
&& rm -f $NAME.raw.bam

# check chromosome content
$SAM_APP idxstats $NAME.bam > $NAME.idxstats.txt

# check duplicate level
echo "=== Check duplication levels ==="
$DEDUP_APP \
--cpu $NCPU \
--in $NAME.bam \
--optical \
--distance $OPTDISTANCE \
--chrskip "chrM|MT"

# check insert size with macs2
echo "=== Macs2 predictd shift prediction ==="
$MACS_APP predictd -i $NAME.bam -g $SIZE --rfile $NAME.predictd


# raw coverage and size prediction with bam2wig
echo "=== bam2wig shift prediction ==="
$BAMWIG_APP \
--in $NAME.bam \
--cpu $NCPU \
--shift \
--model \
--chrom 20 \
--out $NAME



