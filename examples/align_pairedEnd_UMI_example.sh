#!/bin/bash

# for HCI pysano purposes only - delete this and preceding lines and rename as cmd.txt
#e your.name@hci.utah.edu
#c kingspeak
#t 12

# paired-end ChIPSeq with third Fastq file of UMI codes
# alignment with Novoalign followed by duplication and insertion size metrics

# set the base name of the output
NAME=MYNAME

# example mouse Novoalign index
INDEX=/tomato/dev/data/Mouse/Mm10/mm10.nov.illumina.nix
BLACKLIST=/tomato/dev/data/Mouse/Mm10/mm10.blacklist.bed

# optical distance
OPTDISTANCE=2500

## Select the appropriate adaptors

# adapters for standard compatible Illumina TruSeq
ADAPTF=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
ADATPR=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT


# application paths
APP=/tomato/dev/app
NOVO_APP=$APP/novoalign/4.02.02/novoalign
SAM_APP=$APP/samtools/1.16/samtools
MERGE_APP=$APP/modulesoftware/merge_umi_fastq
DEDUP_APP=$APP/modulesoftware/bam_umi_dedup
PICARD_APP=$APP/picard/2.23.3/picard.jar

# merge UMI fastq
echo "=== Merging UMI fastq"
$MERGE_APP *.fastq.gz

# align UMI Fastq
echo "=== Aligning ==="
$NOVO_APP \
-d $INDEX \
-a $ADAPTF $ADATPR \
-C \
-o SAM \
--tune NOVASEQ \
-f *.umi.fastq.gz \
| $SAM_APP fixmate -r -m - $NAME.raw.bam


# sort and index
echo "=== Sorting"
$SAM_APP sort -m 4G -@ 12 -o $NAME.sort.bam $NAME.raw.bam \
&& $SAM_APP index -@ 12 $NAME.sort.bam \
&& rm -f $NAME.raw.bam *.umi.fastq.gz


# remove UMI duplicates
echo "=== Removing UMI duplicates"
$DEDUP_APP \
--in $NAME.sort.bam \
--out $NAME.bam \
--distance $OPTDISTANCE \
--cpu 12 \
&& rm $NAME.sort.bam*


# check chromosome content
$SAM_APP index -@ 12 $NAME.bam
$SAM_APP idxstats $NAME.bam > $NAME.idxstats.txt


# check insertion sizes
echo "=== Check insertion sizes ==="
java -jar $PICARD_APP CollectInsertSizeMetrics \
-I $NAME.bam \
-O $NAME.insertSizeMetrics.txt \
-H $NAME.insertSizes.pdf \
-M 0.5 --VALIDATION_STRINGENCY SILENT --VERBOSITY WARNING 




