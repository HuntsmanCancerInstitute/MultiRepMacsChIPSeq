#!/bin/bash

# run separate paired-end and single-end ChIP calls
# then combine data and do a joint merge peak call

# module load multirepchipseq
export PATH=/Users/tim/test/multirepchipseq_test:$PATH

rm -rf tup1 rpd3 Rpd3_Tup1_merge

echo;echo "====================== Rpd3 only paired-end ======================";echo;
multirep_macs2_pipeline.pl \
--chip Rpd3_Ch1.bam \
--name Rpd3_Ch1 \
--chip Rpd3_Ch2.bam \
--name Rpd3_Ch2 \
--chip Rpd3_Ch3.bam \
--name Rpd3_Ch3 \
--control Rpd3_Input.bam \
--pe \
--dir rpd3 \
--out rpd3_all \
--genome 230000 \
--dupfrac 0.05 \
--optdist 2500 \
--size 275 \
--cutoff 3 \
--peaksize 200 \
--peakgap 100 \
--plot \
--plotpeak /Users/tim/Documents/Code/MultiRepMacsChIPSeq/bin/plot_peak_figures.R \
--cpu 1 \
--job 4 


echo;echo "====================== Tup1 only single-end ======================";echo;
multirep_macs2_pipeline.pl \
--chip Tup1_Ch1.bam \
--name Tup1_Ch1 \
--chip Tup1_Ch2.bam \
--name Tup1_Ch2 \
--chip Tup1_Ch3.bam \
--name Tup1_Ch3 \
--control Tup1_Input.bam \
--dir tup1 \
--out tup1_all \
--genome 230000 \
--dupfrac 0.05 \
--optdist 2500 \
--size 275 \
--cutoff 3 \
--peaksize 200 \
--peakgap 100 \
--plot \
--plotpeak /Users/tim/Documents/Code/MultiRepMacsChIPSeq/bin/plot_peak_figures.R \
--cpu 1 \
--job 4 


echo;echo "====================== Merging Rpd3 and Tup1 Results ======================";echo;

mkdir -p Rpd3_Tup1_merge/Analysis
mkdir Rpd3_Tup1_merge/Count
mkdir Rpd3_Tup1_merge/Fragment
mkdir Rpd3_Tup1_merge/Log2FE
mkdir Rpd3_Tup1_merge/QValue
mkdir Rpd3_Tup1_merge/Peaks

# merge files
join_data_file.pl -o Rpd3_Tup1_merge/Analysis/joint_samples.txt \
rpd3/Analysis/rpd3_all_samples.txt \
tup1/Analysis/tup1_all_samples.txt

join_data_file.pl -o Rpd3_Tup1_merge/Analysis/joint.dedup-stats.txt \
rpd3/Analysis/rpd3_all.dedup-stats.txt \
tup1/Analysis/tup1_all.dedup-stats.txt


# Move bigWigs
mv rpd3/Count/*.bw tup1/Count/*.bw Rpd3_Tup1_merge/Count/
mv rpd3/Log2FE/*.bw tup1/Log2FE/*.bw Rpd3_Tup1_merge/Log2FE/
mv rpd3/Fragment/*.bw tup1/Fragment/*.bw Rpd3_Tup1_merge/Fragment/
mv rpd3/QValue/*.bw tup1/QValue/*.bw Rpd3_Tup1_merge/QValue/


# Run recall
recall_peaks.pl \
--dir Rpd3_Tup1_merge \
--in joint \
--out merge \
--cutoff 3 \
--peaksize 200 \
--peakgap 100 \
--cpu 1 \
--job 4 \
--plot \
--plotpeak /Users/tim/Documents/Code/MultiRepMacsChIPSeq/bin/plot_peak_figures.R 



