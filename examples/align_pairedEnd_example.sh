#!/bin/bash

# for HCI pysano purposes only - delete this and preceding lines and rename as cmd.txt
#e your.name@hci.utah.edu
#c kingspeak
#t 12

# paired-end ChIPSeq, including paired-end ATACSeq
# alignment with Novoalign followed by duplication and insertion size metrics

# set the base name of the output
NAME=MYNAME

# example mouse Novoalign index
INDEX=/tomato/dev/data/Mouse/Mm10/mm10.nov.illumina.nix
BLACKLIST=/tomato/dev/data/Mouse/Mm10/mm10.blacklist.bed

# optical distance
# use 100 for HiSeq runs or 2500 for NovaSeq runs
OPTDISTANCE=2500

## Select the appropriate adaptors

# adapters for standard compatible Illumina TruSeq
ADAPTF=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
ADATPR=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# # adapters for Tn5 ATACSeq
# ADAPTF=CTGTCTCTTATACACATCT
# ADATPR=CTGTCTCTTATACACATCT

# application paths
APP=/tomato/dev/app
NOVO_APP=$APP/novoalign/4.02.02/novoalign
SAM_APP=$APP/samtools/1.16/samtools
DEDUP_APP=$APP/modulesoftware/bam_partial_dedup
PICARD_APP=$APP/picard/2.23.3/picard.jar

# align
echo "=== Aligning ==="
$NOVO_APP \
-d $INDEX \
-a $ADAPTF $ADATPR \
-o SAM \
--tune NOVASEQ \
-f *.gz \
| $SAM_APP fixmate -r -m - $NAME.raw.bam

# sort and index
echo "=== Sorting ==="
$SAM_APP sort -m 4G -@ 8 -o $NAME.bam $NAME.raw.bam \
&& $SAM_APP index -@ 8 $NAME.bam \
&& rm -f $NAME.raw.bam

# check chromosome content
$SAM_APP idxstats $NAME.bam > $NAME.idxstats.txt

# check duplicate level
echo "=== Check duplication levels ==="
$DEDUP_APP --pe --cpu 12 \
--in $NAME.bam \
--optical --distance $OPTDISTANCE \
--blacklist $BLACKLIST \
--chrskip 'chrM|MT' \
--report


# check insertion sizes
echo "=== Check insertion sizes ==="
java -jar $PICARD_APP CollectInsertSizeMetrics \
-I $NAME.bam \
-O $NAME.insertSizeMetrics.txt \
-H $NAME.insertSizes.pdf \
-M 0.5 --VALIDATION_STRINGENCY SILENT --VERBOSITY WARNING 




