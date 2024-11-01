# Transcription Factor ChIPSeq Partial Duplication Example

An example of a DNA-binding transcription factor ChIPSeq from mammalian cells. The
ChIP samples were prepared with NEBNext Ultra II DNA Library kit and sequenced on a
NovaSeq 6000. This ChIPSeq had a total of 108M alignments. After removing 15M optical
duplicates (based on a pixel distance of 2500), the UMI duplication rate was 36%
resulting in 60M unique alignments. This was one of three comparable biological
duplicates.

See the [BASH script](https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq/blob/master/docs/DeDuplicationEvaluation/duplication_comparison_cmd.sh)
for the exact pipeline used in processing.

### Peak call number

This transcription factor has relatively few occupied sites, calling only about 1500 
sites in the genome. There is very little difference between the number of peaks 
called with UMI-deduplicated and full de-duplication, and adding in a subsample of 
duplicate alignments has a modest effect, reducing or increasing by a small number. 
Leaving in all duplicates calls nearly 40% more peaks. 

The mean peak length was 418 ± 182 bp, median 370 bp at a threshold q-value of 2. 

![TF_peak_number](DeDuplicationEvaluation/TF.peak_number.png)


### ChIP Efficiency

The fraction of alignments within the respectively called peaks was quite low, as 
expected for such few peak numbers called, ranging from 0.25% to 0.30%.  

![TF_chip_efficiency](DeDuplicationEvaluation/TF.chip_efficiency.png)


### Peak Fragment Coverage Profile

The depth-normalized fragment coverage profile over the called peak midpoint ± 1 Kb 
shows no difference between the subsets.

![TF_profile_fragment](DeDuplicationEvaluation/TF_profile_fragment_hm.png)


### Comparison of Duplicate Alignments Within Peaks

To compare the influence of retained duplicate alignments on the ChIP efficiency, the
fraction of alignments (used in the peak call) within each respectively called peaks
was compared with the fraction of completely de-duplicated alignments. Each fraction
was plotted below.

In all cases, a higher percentage of duplicate-containing alignments were observed in
the peaks than completely deduplicated alignments, indicating there are indeed
duplicate alignments found within peaks. Notably, however, the difference is quite
small. 

![TF_efficiency_comparison](DeDuplicationEvaluation/TF_comparison.chip_efficiency.png)


### Conclusion

In this transcription factor ChIPSeq, peak calling without any duplicate alignments
most closely matches the true peak call set with UMI deduplicated alignments when
considering the number of peaks called. However, both 5% and 10% duplication levels
generated very similar numbers, indicating minimal effect of duplication alignments
on call sets.

When comparing the fraction of alignments with and without duplicates within called
peaks, there is very little additional fraction of duplicates within peaks. In fact,
the percentages are almost identical with UMI de-duplication, indicating that
biological duplication does not in fact play a large contributing factor.




