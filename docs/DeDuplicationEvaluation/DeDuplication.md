# MultiRepMacsChIPSeq - De-Duplication

|[Home](../Readme.md)|[Overview](../Overview.md)|[Usage](../Usage.md)|[Variations](../Variations.md)|[Examples](../Examples.md)|[Applications](../applications.md)|[Install](../Install.md)|

## De-Duplication sub-sampling evaluation

The Multi-Replica Macs ChIPSeq Wrapper includes a novel application,
[bam_partial_dedup](applications/bam_partial_dedup.md), for sub-sampling
duplicate alignments to a standard fractional level across all replicates and
samples. The theory behind the application is to keep some fraction of duplicate
alignments on the assumption that these may include biological duplicates, known to
exist at peaks of strong enrichment, along with PCR and technical duplication. By
normalizing the fraction of duplicates, as well as removing known sources of
technical duplication, including optical duplicates and empirically-derived high-copy
number genomic intervals, peak calling, and more importantly differential analysis,
would be improved over standard methods of peak calling (where all duplicates are
excluded).

Test of this hypothesis was hampered by knowing exactly what was biological versus
amplification artifact. Optical duplicates could generally be classified by their
physical distance from their nearest duplicate neighbor on an Illumina flow cell, but
was still arbitrary at best. The use of Unique Molecular Indexes (UMIs) changes this.
By including a random nucleotide sequence within the adapter used in library
preparation, each molecule can be uniquely tagged prior to amplification. This tag,
in conjunction with alignment coordinates, can differentiate between an amplification
or biological duplicate simply by comparing the UMI sequence.

Here, we test the hypothesis that random duplicate subsampling can improve signal and
peak calls over complete de-duplication or no de-duplication. ChIPSeq analysis was
performed against histone mark H3K4Me3, RNA Polymerase II, and a bHLH Transcription
Factor. 


### Methods

Libraries were prepared using NEBNext Ultra II DNA library kit, which includes a
random UMI sequence in the adapter. Libraries were sequenced on an Illumina Novaseq
Novaseq 6000. The UMI was sequenced as third index read (11 bases), and recorded as
a thirdFastq file during de-multiplexing with Illumina Bcl2Fastq software.

The UMI fastq was appended to the two primary reads as SAM attribute tag `RX` using
the
[merge_umi_fastq.pl](https://huntsmancancerinstitute.github.io/UMIScripts/apps/merge_umi_fastq.html)
application from the
[UMIScripts package](https://huntsmancancerinstitute.github.io/UMIScripts).
Reads were aligned to the genome using Novocraft Novoalign (version 4.3).

Alignment bam files were processed with either
[bam_partial_dedup.pl](../applications/bam_partial_dedup) or with
[bam_umi_dedup.pl](https://huntsmancancerinstitute.github.io/UMIScripts/apps/bam_umi_dedup.html)
from the [UMIScripts package](https://huntsmancancerinstitute.github.io/UMIScripts/).
With partial de-duplication, all optical duplicates (optical threshold of 2500
pixels) were completely removed. Exclusion peaks, identified from the corresponding
Input, were called (without prior de-duplication) and used to exclude alignments in
all cases.

Six different permutations were generated and compared. Labels are with respect to
the number of duplicates retained.

- No duplicates retained (label `no`)
- Duplicates subsampled to 5 percent (label `p05`)
- Duplicates subsampled to 10 percent (label `p10`)
- Duplicates subsampled to 20 percent (label `p20`)
- True duplicates removed using UMIs (label `umi`)
- All duplication retained (label `full`)

The processed bam files were then put through the MultiRepChIPSeq pipeline (version
17.9) to call peaks as separate samples, merged, and re-scored. The complete pipeline
(for those wanting to follow along at home) is found in the adjacent [BASH
script](https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq/blob/master/docs/DeDuplicationEvaluation/duplication_comparison_cmd.sh). 

Sequence files are from unpublished results and currently unavailable.


### Results

Several plots are provided here as generated by the pipeline.

See [H3K4Me3 results](H3K4Me3_results.md) for a high
occurrence, high occupancy factor.

See [RNA Polymerase II results](RNAPol2_results.md) for a
very high occurrence, moderate occupancy factor.

See [Transcription Factor results](TF_results.md) for a low
occurrence, low occupancy factor.


### Conclusion

Including a sub-sampling of duplicate alignments does not necessarily replicate the
results obtained with UMI-specific de-duplication. In some cases, UMI de-duplication
is closer to full de-duplication (no replicates), while in others it is closer to low
partial de-duplication (5 percent). Increasing the amount of duplicates, up to the
full amount (including optical duplicates), results in an increasing number of peaks
called, presumably false positive.

There is a modest increase in the counts and fraction of retained duplicates within
called peaks. However, due to depth-normalization of the samples, this does not
result in significant change in fragment pileup or log2 fold enrichment.

In most cases, including a partial fraction of alignment duplicates has a negligible
to modest effect. The broader the expected peak, the less benefit is observed. For
broad domains, such as H3K27 modifications, no benefit is expected (data not shown).

These results are applicable for sonicated (random fragmentation) chromatin. For
applications with enzymatic fragmentation, including MNase and ATACSeq, these results
may not be applicable. Enzymatic fragmentation frequently results in considerably
higher rates of biological duplication due to the restrictive nature of enzyme
accessibility to genomic DNA. Rates of duplication in excess of 50% is common, not
counting optical duplicates. Small genomes, such as yeast, exacerbate this even
further with rates upwards of 90%. 

Without a doubt, using UMI-codes should be the gold standard for handling duplicates.
Until this is standard, using workarounds, such as duplicate sub-sampling, has
marginal benefit in most cases, depending on the scenario. Keeping all
duplicates, especially with high optical duplicate levels (usually generated by
machines with patterned flowcells such as Illumina Nextseq or Novaseq), is
detrimental and can lead to false positives. Discarding all duplicates, the
traditional standard, is still ideal when fragments are randomly generated, e.g.
shearing by sonication, but caution should be considered with experiments using
enzymatically digested non-random fragmentation, e.g. ChIP of MNase fragments. 

Hence, factors to consider include:

1. level of optical duplicates (which sequencing machine was used)
2. method of fragmentation
3. enrichment pattern: narrow or broad

Note that sequencing exported from SRA do not always keep original read names and
optical duplicates cannot be distinguished. Therefore, full de-duplication is the only 
option with SRA exported reads.


