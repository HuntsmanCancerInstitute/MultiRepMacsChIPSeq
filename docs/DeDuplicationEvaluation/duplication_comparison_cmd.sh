#!/bin/bash

# script to run the multi-sample ChIP pipeline to test different
# levels of de-duplication, including with UMIs.

# software packages are installed under Gnu modules
module load multirepchipseq umiscripts parallel

# substitute MYCHIP and MYINPUT with appropriate sample names
CHIP=MYCHIP
INPUT=MYINPUT


echo "========== UMI de-duplication =========="
echo "===== ChIP"
bam_umi_dedup.pl \
--in $CHIP.bam \
--out umi_ChIP.bam \
--distance 2500 --cpu 12 \
--umi RX > umi_ChIP.out.txt

echo "===== Input"
bam_umi_dedup.pl \
--in $INPUT.bam \
--out umi_In.bam \
--distance 2500 --cpu 12 \
--umi RX > umi_In.out.txt


echo;echo "========== Complete de-duplication =========="
echo "===== ChIP"
bam_partial_dedup.pl \
--in $CHIP.bam \
--out no_ChIP.bam \
--distance 2500 --cpu 12 --chrskip 'chrM' --pe \
--blacklist exclusion.bed \
--frac 0 --max 1 > no_ChIP.out.txt

echo "===== Input"
bam_partial_dedup.pl \
--in $INPUT.bam \
--out no_In.bam \
--distance 2500 --cpu 12 --chrskip 'chrM' --pe \
--blacklist exclusion.bed \
--frac 0 --max 1 > no_In.out.txt

echo;echo "========== Partial de-duplication at 5 percent =========="
echo "===== ChIP"
bam_partial_dedup.pl \
--in $CHIP.bam \
--out 05p_ChIP.bam \
--distance 2500 --cpu 12 --chrskip 'chrM' --pe \
--blacklist exclusion.bed \
--frac 0.05 --seed 1 > 05p_ChIP.out.txt

echo "===== Input"
bam_partial_dedup.pl \
--in $INPUT.bam \
--out 05p_In.bam \
--distance 2500 --cpu 12 --chrskip 'chrM' --pe \
--blacklist exclusion.bed \
--frac 0.05 --seed 1 > 05p_In.out.txt

echo;echo "========== Partial de-duplication at 10 percent =========="
echo "===== ChIP"
bam_partial_dedup.pl \
--in $CHIP.bam \
--out 10p_ChIP.bam \
--distance 2500 --cpu 12 --chrskip 'chrM' --pe \
--blacklist exclusion.bed \
--frac 0.1 --seed 1 > 10p_ChIP.out.txt

echo "===== Input"
bam_partial_dedup.pl \
--in $INPUT.bam \
--out 10p_In.bam \
--distance 2500 --cpu 12 --chrskip 'chrM' --pe \
--blacklist exclusion.bed \
--frac 0.1 --seed 1 > 10p_In.out.txt

echo;echo "========== Partial de-duplication at 20 percent =========="
echo "===== ChIP"
bam_partial_dedup.pl \
--in $CHIP.bam \
--out 20p_ChIP.bam \
--distance 2500 --cpu 12 --chrskip 'chrM' --pe \
--blacklist exclusion.bed \
--frac 0.2 --seed 1 > 20p_ChIP.out.txt

echo "===== Input"
bam_partial_dedup.pl \
--in $INPUT.bam \
--out 20p_In.bam \
--distance 2500 --cpu 12 --chrskip 'chrM' --pe \
--blacklist exclusion.bed \
--frac 0.2 --seed 1 > 20p_In.out.txt


echo;echo "========== Combining duplication stats"
combine_std_chipstats.pl duplication_stats.txt *.out.txt

parallel -k -j 1 echo "======= {} =======" '>>' deduplication.logs.txt \
 '&&' cat {} '>>' deduplication.logs.txt \
 ':::' *.out.txt
rm *.out.txt


# link original bam files to comply with naming scheme
ln -s $CHIP.bam full_ChIP.bam
ln -s $CHIP.bam.bai full_ChIP.bam.bai
ln -s $INPUT.bam full_In.bam
ln -s $INPUT.bam.bai full_In.bam.bai


echo;echo "===================== ChIP pipeline of all at q-value 2 ====================="
multirep_macs2_pipeline.pl \
--chip no_ChIP.bam \
--control no_In.bam \
--name no \
--chip 05p_ChIP.bam \
--control 05p_In.bam \
--name 05p \
--chip 10p_ChIP.bam \
--control 10p_In.bam \
--name 10p \
--chip 20p_ChIP.bam \
--control 20p_In.bam \
--name 20p \
--chip full_ChIP.bam \
--control full_In.bam \
--name full \
--chip umi_ChIP.bam \
--control umi_In.bam \
--name umi \
--dir dup_comparison \
--out all_q2 \
--pe \
--nodedup \
--blacklist exclusion.bed \
--cutoff 2 \
--peaksize 300 \
--peakgap 100 \
--chrskip 'chrM' \
--cpu 8 \
--job 4

echo;echo "===================== Recall peaks at q-value 3 ====================="
recall_peaks.pl \
--dir dup_comparison \
--in all_q2 \
--out all_q3 \
--cutoff 3 \
--peaksize 300 \
--peakgap 100 \
--cpu 8 \
--job 4



echo;echo "===================== Marking remaining duplicates ====================="

parallel -j 2 -v \
bam_partial_dedup.pl \
--in {1}_{2}.bam \
--out {1}_{2}.mark.bam \
--blacklist exclusion.bed \
--pe --max 1 --mark --cpu 8 ">" {1}_{2}.mark.out.txt \
':::' 05p 10p 20p umi full \
':::' ChIP In

combine_std_chipstats.pl used_duplication_stats.txt *.mark.out.txt

parallel -k -j 1 echo "======= {} =======" '>>' used_duplication_deduplication.logs.txt \
 '&&' cat {} '>>' used_duplication_deduplication.logs.txt \
 ':::' *.out.txt
rm *.out.txt


echo;echo "===================== Extracting remaining duplicates ====================="

SAM=`which samtools`

echo; echo "==== indexing"
parallel -j 2 -v \
$SAM index -@ 8 {1}_{2}.mark.bam \
':::' 05p 10p 20p umi full \
':::' ChIP In

echo; echo "==== extracting"
parallel -j 2 -v \
$SAM view --threads 8 -f DUP,PROPER_PAIR -o {1}_{2}.dup.bam {1}_{2}.mark.bam \
':::' 05p 10p 20p umi full \
':::' ChIP In

echo; echo "==== indexing"
parallel -j 2 -v \
$SAM index -@ 8 {1}_{2}.dup.bam \
':::' 05p 10p 20p umi full \
':::' ChIP In


echo;echo "===================== Generating Duplicate Count data ====================="

mkdir dup_comparison/DupCount dup_comparison/DupEfficiency

parallel -j 2 -v \
bam2wig.pl \
--bw --pe --mid --cpu 8 \
--blacklist exclusion.bed \
--chrskip chrM  \
--in {1}_{2}.dup.bam \
--out dup_comparison/DupCount/{1}_{2}.count.bw \
':::' 05p 10p 20p umi full \
':::' ChIP In




echo;echo "===================== Collecting Duplicate Efficiency ====================="

parallel -j 2 -v \
get_chip_efficiency.pl \
--in dup_comparison/Peaks/{}.bed \
--group dup_comparison/Analysis/all_q2_samples.txt \
--out dup_comparison/{}.dup_efficiency.txt  \
dup_comparison/DupCount/{}_ChIP.count.bw \
dup_comparison/DupCount/{}_In.count.bw \
':::' 05p 10p 20p full umi

join_data_file.pl \
--out dup_comparison/DupEfficiency/all_q2.chip_efficiency.txt \
dup_comparison/*.dup_efficiency.txt

rm dup_comparison/*.dup_efficiency.txt

parallel -j 2 -v \
get_chip_efficiency.pl \
--in dup_comparison/Peaks1/{}.bed \
--group dup_comparison/Analysis1/all_q3_samples.txt \
--out dup_comparison/{}.dup_efficiency.txt  \
dup_comparison/DupCount/{}_ChIP.count.bw \
dup_comparison/DupCount/{}_In.count.bw \
':::' 05p 10p 20p full umi

join_data_file.pl \
--out dup_comparison/DupEfficiency/all_q3.chip_efficiency.txt \
dup_comparison/*.dup_efficiency.txt

rm dup_comparison/*.dup_efficiency.txt


echo;echo "===================== Plotting Duplicate Efficiency ====================="

plot_peak_figures.R -i dup_comparison/DupEfficiency/all_q2
plot_peak_figures.R -i dup_comparison/DupEfficiency/all_q3



# Cleanup
rm *.mark.bam*
rm full_ChIP.bam* full_In.bam*




