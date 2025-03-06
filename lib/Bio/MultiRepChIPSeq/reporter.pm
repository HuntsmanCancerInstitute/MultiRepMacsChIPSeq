package Bio::MultiRepChIPSeq::reporter;

use strict;
use English qw(-no_match_vars);
use Carp;
use File::Spec::Functions qw( catfile splitpath );
use Bio::ToolBox 1.70;
use Bio::ToolBox::utility qw(format_with_commas);

our $VERSION = 20.1;

sub add_header_report {
	my $self = shift;
	my $provided_options = shift;

	my $title = sprintf "%s %s", $self->dir, $self->out;
	my $v     = $self->version;
	my $log   = $self->out . "_job_output_logs.txt";
		
	# header
	my $string = <<END;
# Pipeline for $title

This is a report from the 
[Multi-Replica Multi-Sample MACS2 ChIPSeq pipeline](https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq),
Version $v.

Individual application outputs may be found in the combined log output file `$log`.
END

	# folder descriptions
	if ($self->organize) {
		my $suffix = $self->dir_suffix;
		$string .= <<END;

The output is organized into the following folders


- **Analysis$suffix**: Collected data and QC files. All files are text, with data files
  gzip compressed.
END
		if ( $self->{progress}{fragment} ) {
			$string .= <<END;
- **Count**: BigWig files representing depth-normalized, scaled, fragment counts at 
  the fragment mid point (or cut-site end point) coordinate, i.e. single-point data.
- **Fragment**: BigWig files representing depth-normalized (per Million) fragment
  coverage.
- **Log2FE**: BigWig files representing Log 2 scaled Fold Enrichment (Signal over
  lambda control reference or Input).
END
		}
		$string .= <<END;
- **Peaks$suffix**: The called-peak, genomic-interval files in either narrowPeak,
  gappedPeak (broad calls only), or bed (6-column, merged peaks only) file formats.
- **PeakSummits$suffix**: Bed files representing the peak signal summit coordinate
  (1 bp) or peak midpoint (merged peaks only) for each peak call.
END
		if ($self->independent) {
			$string .= <<END;
- **Replicate_Plots**: Intesection and QC Plots for independently called replicate peaks.
- **Replicate-Mean_Plots**: Intersection, Analysis, and QC plots for the called peaks
  from the depth-normalized, replicate-average fragment coverage.
- **Replicate-Merge_Plots**: Intersection, Analysis, and QC plots for the merged,
  independently called replicate peaks.
END
		}
		else {
			$string .= <<END;
- **Plots$suffix**: Intersection, Analysis, and QC plots for the called sample peaks.
END
		}
		if ( $self->{progress}{fragment} ) {
			$string .= <<END;
- **QValue**: BigWig files representing the multiple-test corrected statistical
  enrichment (poission distribution) over each base pair in the genome. Values are
  transformed as `-1*log10` for ease in displaying genome browsers.
END
		}
		if ( $self->savebam and $self->dedup ) {
			$string .= <<END;
- **DeDupBam**: Saved de-duplicated bam files used in the analysis.
END
		}
		if ( $self->savebdg ) {
			$string .= <<END;
- **BedGraph**: Saved intermediate bedGraph files prior to bigWig conversion.
END
		}
	}

	# pipeline parameters section
	$string .= <<END;
### Parameters

The pipeline was run with the following explicit parameters:

	$PROGRAM_NAME $provided_options

END

	return $string;
}

sub add_samples_report {
	my $self = shift;

	my $string = <<END;
### Samples

The following samples (conditions) and replicate files were provided as input to
pipeline.

END

	# samples
	my $i = 1;
	foreach my $Job ( $self->list_jobs ) {
		$string .= sprintf( "%d. %s ChIP: %s\n", $i, $Job->job_name,
			join( ', ', map { sprintf "`%s`", $_ } ( $Job->chip_bams ) ) );
		$i++;
		$string .= sprintf( "%d. %s Control: %s\n", $i, $Job->job_name,
			join( ', ', map { sprintf "`%s`", $_ } ( $Job->control_bams ) ) );
		$i++;
	}

	# universal control
	if ( $self->has_universal_control ) {
		$string .= <<END;

Only one reference control was provided, which was used as a universal control
for all provided samples.
END
	}

	$string .= <<END;

--------

## Pipeline Options

END

	return $string;
}

sub add_genome_report {
	my $self = shift;
	my $string;
	
	# genome
	if ( $self->genome and exists $self->{all_map_fraction} ) {
		my $size      = format_with_commas( $self->genome );
		my $full      = format_with_commas( $self->{full_genome_size} );
		my $all_frac  = sprintf "%s%%", $self->{all_map_fraction};
		my $uniq      = format_with_commas( $self->{unique_map} );
		my $uniq_frac = sprintf "%s%%", $self->{unique_map_fraction};
		my $exclusion = $self->chrskip;
		$string .= <<END;

### Genome Size

The full genome size, minus excluded chromosomes (`$exclusion`), is $full bp.

All provided alignments were empirically determined to cover $size bp,
or `$all_frac`. Uniquely mapping alignments cover $uniq bp, or $uniq_frac.

For enrichment calculations, $size bp was used.

END
		if ( $all_frac < 75 ) {
			$string .= <<END;
**WARNING**: The mapped fraction of the genome is low and may potentially
negatively impact background estimation and q-value scores. Manually setting
the genome size may improve peak calling.

END
		}
	}
	elsif ( $self->genome ) {
		my $size = format_with_commas( $self->genome );
		my $full = format_with_commas( $self->{full_genome_size} );
		my $fraction = sprintf "%.1f%%",
			( $self->genome / $self->{full_genome_size} ) * 100;
		my $exclusion = $self->chrskip;
		$string .= <<END;

### Genome Size

The full genome size, minus excluded chromosomes (`$exclusion`), is $full bp.

The provided genome size of $size bp, or $fraction, was used for enrichment
calculations. 

END
	}

	return $string;
}

sub add_exclusion_report {
	my $self = shift;
	my $string;

	# exclusion list
	if ( $self->{progress}{control_peak} ) {
		# input generated exclusion list
		my $file = ( splitpath( $self->exclude ) )[2];
		my ($control, $count);
		if ( $self->organize ) {
			$control  = catfile( 'Peaks', $file );
			$count = format_with_commas( $self->count_file_lines( catfile( $self->dir,
				$control ) ) );
		}
		else {
			$control  = $file;
			$count = format_with_commas( $self->count_file_lines( $self->exclude ) );
		}
		$string .= <<END;

### Exclusion List

An exclusion list (also known as gray list) was generated from the provided
reference control (Input) Bam files. Peaks were called from the combined
coverage of all control Bam files using MACS2. These typically represent
repetitive sequences (simple repeats or repetitive sequences), often near
centromeres, telomeres, or other heterochromatic regions. Alignments overlapping
these intervals will be excluded from futher analysis. This reduces a source of
duplicate and non-informative alignments and minimizes false positive peak
calls.

There were **$count** exclusion intervals identified in the file `$control`.
END
		
	}
	elsif ( $self->exclude ) {
		# user supplied
		my $f = $self->exclude;
		$string .= <<END;

### Exclusion List

The user-provided exclusion interval file, `$f`, was used to exclude alignments
from the Bam files.
END
	}

	return $string;
}

sub add_deduplication_report {
	my $self = shift;

	my $fraction = $self->dupfrac;
	my $maxdepth = $self->maxdepth || 0;
	my $optdist  = $self->optdist;
	my $pair     = $self->deduppair;
	
	my $string   = <<END;

### Alignment De-Duplication

Prior to analysis, duplicate alignments were removed from provided bam files.
END
		
	# optical distance
	if ( $optdist > 0 ) {
		$string .= <<END;
Optical duplicates (sequencing artifacts) were identified with a pixel distance 
less than $optdist and discarded.
END
	}

	# method of deduplication
	if ( $fraction > 0 and $maxdepth > 0 ) {
		$string .= <<END;
Non-optical duplicates were sub-sampled until the duplication rate was $fraction
or less and the maximum depth at any position was $maxdepth or less.
END
	}
	elsif ( $fraction > 0 and $maxdepth == 0 ) {
		$string .= <<END;
Non-optical duplicates were sub-sampled until the duplication rate was $fraction
or less.
END
	}
	elsif ( $fraction == 0 and $maxdepth > 0 ) {
		$string .= <<END;
Duplicate alignments were removed until the maximum depth at any position was 
$maxdepth or less.
END
	}

	# alignment count summary
	my $data = catfile( 'Analysis', $self->out . '.dedup-stats.txt' );
	$string .= <<END;

Alignment counts are in the file `$data`.
END
	if ($self->plot) {
		my $plot;
		if ($self->independent) {
			$plot = catfile( 'Replicate-Mean_Plots', $self->out .
				'.duplicate-counts.png' );
		}
		else {
			$plot = catfile( 'Plots', $self->out . '.duplicate-counts.png' );
		}
		$string .= <<END
A summary of retained and discarded alignments for each bam file is shown below.

![duplicate_counts]($plot)
END
	}

	return $string;
}

sub add_filter_report {
	my $self = shift;

	# bam filter
	my $string = <<END;

### Alignment Filtering

Unwanted alignments were filtered from the provided bam files. These included marked
duplicates, secondary, supplemental, QC fail, and unmapped alignments.
END
	if ($self->paired) {
		$string .= "Alignment pairs with one unmapped end were discarded.\n";
	}
	if ($self->mapq) {
		my $m = $self->mapq;
		$string .=
			"Alignments with a mapping quality score less than $m were discarded\n";
	}

	return $string;
}

sub add_peak_call_parameters {
	my $self = shift;
	
	my $peaksize = $self->peaksize;
	my $peakgap  = $self->peakgap;
	my $cutoff   = $self->cutoff;
	my $broadcut = $self->broadcut;
	my $broadgap = $self->broadgap;
	my ($decimal, $decimal2);
	{
		my $number  = int($cutoff) + 1;
		my $pattern = qq(%.$number) . 'f';
		$decimal    = sprintf $pattern, ( 1 / ( 10 ** $cutoff ) );
		$number     = int($broadcut) + 1;
		$pattern    = "%.$number" . 'f';
		$decimal2   = sprintf $pattern, ( 1 / ( 10 ** $broadcut ) );
	}
	my $string = <<END;

### Peak Calling parameters

Narrow peaks calls were made using a q-value threshold of **$cutoff** (essentially
a False Discovery Rate <= $decimal), a minimum peak size of **$peaksize** bp, and a
maximum allowed gap of **$peakgap** bp. 
END

	if ($self->broad) {
		$string .= <<END;

For broad, or gapped-peak, calling, significant peaks were called using a q-value
threshold of **$cutoff** (essentially FDR <= $decimal) and maximum gap of
**$peakgap** bp. Dispersed or nearby peaks were merged with a cutoff of
**$broadcut** (essentially FDR <= $decimal2). The maximum linking between
significant peaks was **$broadgap** bp. The minimum peak length was
**$peaksize** bp.
END
	}

	if ($self->independent) {
		$string .= <<END;

Two sets of peak calls were generated using two separate methods, which are
comparable but mutually exclusive. 

- Replicate-Merge
	
	Peak calls are identified separately for each individual replicate and then 
	merged into a single call set for each sample condition. This tolerates
	individual replicate efficiency but requires consensus before merging,
	avoiding sporadic calls made from individual replicates.

- Replicate-Mean

	Individual replicates are averaged together in a depth-normalized manner
	and peak calls generated from the mean signal. This is a consensus peak
	call, as each replicate must contribute signal to the mean and thus can
	avoid some outlier effects.

Depending on the strength and efficiency of each replicate and how well they
correlate with each other, these peak call sets may be nearly identical or very
different. A comparison between the two is shown at the end of the report.
END
	}
	else {
		$string .= <<END;

If replicates were included in sample conditions, these were averaged together
in a depth-normalized manner prior to generating peak calls.
END

	}

	return $string;
}

sub add_coverage_report {
	my $self = shift;

	my $status = $self->paired ? 'paired-end' : 'single-end';
	my $map    = $self->mapq;
	my $string = <<END;

--------

## Fragment Coverage

Fragment coverage tracks were generated from $status alignments, ignoring secondary,
supplementary, QC fail, and marked-duplicate alignments, with a mapping quality value
greater than $map. Alignments overlapping exclusion intervals or on excluded
chromosomes were skipped.
END

	# fragment generation based on single or paired end
	if ($self->paired) {
		my $min  = $self->minsize;
		my $max  = $self->maxsize;
		$string .= <<END;
Fragment coverage was generated by using the endpoints from properly paired alignments
with an insertion size between $min and $max bp. 
END
	}
	elsif ( $self->shiftsize and $self->fragsize ) {
		my $shift = $self->shiftsize;
		my $size  = $self->fragsize;
		$string  .= <<END;
Fragment coverage was generated by shifting the alignment start coordinate by $shift
bp and then extending by $size bp. 
END
	}
	else {
		my $size = $self->fragsize;
		$string .= <<END;
Fragment coverage was generated by extending $size bp from the alignment start point
to approximate the mean-size fragment coverage. 
END
	}

	# multiple hit scaling
	if ( $self->fraction ) {
		$string .= <<END;
Multi-mapping fragment coverage was scaled down proportionally based on the reported
number of hits (SAM flag `NH`) in the genome.
END
	}

	# add fragment counts for each bam file
	$string .= <<END;
Coverage tracks were scaled to a normalized depth of one million fragments (or Reads 
Per Million or RPM). Multiple sample replicates were averaged after scaling.

The following number of fragments were accepted for analysis after filtering.

| File | Number (Millions) |
|---|---|
END

	# add table items
	foreach my $Job ($self->list_jobs) {
		foreach my $file ( $Job->chip_use_bams ) {
			$string .= sprintf "| `%s` | %.3f |\n", $file, $self->{bam2count}{$file};
		}
		foreach my $file ( $Job->control_use_bams ) {
			$string .= sprintf "| `%s` | %.3f |\n", $file, $self->{bam2count}{$file};
		}
	}

	# add count scaling factors
	my $method = $self->{original_depth};
	my $target = sprintf "%.3f", $self->targetdepth;
	if ($method =~ /mean | median | min/x) {
		$string .= <<END;

The $method fragment depth was calculated as $target Million. All count data bigWig
tracks were scaled to this depth.
END
	}
	else {
		# user specified target depth, the calculated median depth was inserted
		# into the original_depth key for reporting purposes
		$method  = sprintf "%.3f", $method;
		$string .= <<END;

The median fragment depth was calculated as $method Million. All count data bigWig
tracks were scaled to the provided target depth of $target Million.
END
	}
	
	# add enrichment scaling factors
	if ( $self->samedepth ) {
		$string .= <<END;

The enrichment signal was calculated for each sample using the provided target depth of
$target Million.
END
	}
	else {
		$string .= <<END;

The following replicate-mean fragment depths were used for calculating enrichment for
each sample. The minimum of a ChIP - Control pair was used; see the log file.

| Sample | Item | Mean Depth |
|---|---|---|
END

		# add items to table
		foreach my $Job ( $self->list_jobs ) {
			if ( $Job->chip_bams ) {
				$string .= sprintf "| %s | ChIP | %.3f |\n", $Job->job_name, 
					$self->seq_depth_for_file( $Job->chip_bdg );
			}
			if ( $Job->control_bams ) {
				$string .= sprintf "| %s | Control | %.3f |\n", $Job->job_name, 
					$self->seq_depth_for_file( $Job->lambda_bdg );
			}
		}
	}

	return $string;
}

sub add_independent_peak_calls_report {
	my $self = shift;

	my $string = <<END;

--------

## Independent Replicate Peak Calls

Peaks calls were made independently for each sample replicate as described
above. The following number of peaks were identified. 

| Sample | Replicate | File | Number |
|---|---|---|---|
END
	foreach my $Job ( $self->list_jobs ) {
		next unless $Job->repmerge_peak;
		my @names = $Job->chip_rep_names;
		my @peaks = $Job->rep_peaks;
		for my $i ( 0 .. $#names ) {
			$string .= sprintf "| %s | %s | `%s` | %s |\n",
				$Job->job_name,
				$names[$i],
				( splitpath( $peaks[$i] ) )[2],
				format_with_commas( 
					$self->count_file_lines( $peaks[$i] )
				);
		}
	}

	if ( $self->plot ) {
		$string .= "\nThe UpSet intersection plots for each sample are shown below.\n";
		foreach my $Job ( $self->list_jobs ) {
			next unless scalar( $Job->chip_use_bams ) > 1;
			my $name = $Job->job_name;
			my $plot = catfile( 'Replicate_Plots', $Job->job_name .
				'.rep_merge.intersection_upset.png' );
			$string .= <<END

#### Replicates for $name Peaks

![$name\_replicate_upset]($plot)
END
		}
	}

	return $string;
}

sub add_independent_broad_peak_calls_report {
	my $self = shift;

	my $string = <<END;

--------

## Independent Gapped Replicate Peak Calls

Broad, gapped peaks were called independently for each sample replicate as
described above. The following number of broad peaks were identified.

| Sample | Replicate | File | Number |
|---|---|---|---|
END
	foreach my $Job ( $self->list_jobs ) {
		next unless $Job->repmerge_gappeak;
		my @names = $Job->chip_rep_names;
		my @peaks = $Job->rep_gappeaks;
		for my $i ( 0 .. $#names ) {
			$string .= sprintf "| %s | %s | `%s` | %s |\n",
				$Job->job_name,
				$names[$i],
				( splitpath( $peaks[$i] ) )[2],
				format_with_commas( 
					$self->count_file_lines( $peaks[$i] )
				);
		}
	}

	if ( $self->plot ) {
		$string .= 
			"\nThe UpSet intersection plots for each sample are shown below.\n";
		foreach my $Job ( $self->list_jobs ) {
			next unless scalar( $Job->chip_use_bams ) > 1;
			my $name = $Job->job_name;
			my $plot = catfile( 'Replicate_Plots', $Job->job_name .
				'.rep_merge.intersection_upset.png' );
			$string .= <<END

#### Replicates for $name Broad Peaks

![gap_$name\_replicate_upset]($plot)
END
		}
	}

	return $string;
}


sub add_merged_replicates_report {
	my $self = shift;

	my $overlap = $self->minpeakover || 'n-1';
	my $gap     = $self->fragsize;
	my $rep_num = 0;  # total number of replicates
	my $string  = <<END;

--------

## Merged-replicate Peaks

Independently called replicate peaks were merged for each sample, requiring at least
$overlap overlaps and a maximum gap of $gap bp. The final number of peaks for each
sample are listed below.

| Sample | File | Number |
|---|---|---|
END
	foreach my $Job ($self->list_jobs) {
		next unless $Job->repmerge_peak;
		$string .= sprintf "| %s | `%s` | %s |\n", $Job->job_name,
			( splitpath( $Job->repmerge_peak ) )[2],
			format_with_commas( $self->count_file_lines( $Job->repmerge_peak ) );
		$rep_num += scalar( $Job->chip_rep_names );
	}
	
	# Merged replicates
	my $merged_base  = ( splitpath( $self->repmerge_merge_base ) )[2];
	my $merged_peak  = $merged_base . '.bed';
	my $plot_dir     = 'Replicate-Merge_Plots';
	my $merged_count = format_with_commas( $self->count_file_lines(
		$self->repmerge_merge_base . '.bed' ) );
	$string .= <<END;

These merged-sample peaks were merged into a final master replicate-merge peak set
with **$merged_count** intervals. This file, `$merged_peak`, was then used for sample
comparison.
END

	# add plots
	if ( $self->plot and $merged_count ) {
		my $merge_upset  = catfile( $plot_dir, $merged_base . '.intersection_upset.png' );
		my $merge_length = catfile( $plot_dir, $merged_base . '.peak_lengths.png' );
		$string .= <<END;
The UpSet intersection between samples and sample peak length distributions are shown
below. Additional plots may be found in the **$plot_dir** folder.

![rep-merge_upset]($merge_upset)

![rep-merge_peak_length]($merge_length)
END
	}

	# Merged replicates correlation
	if ( $self->plot and length($merged_count) > 1 and $rep_num > 2 ) {
		my $pca      = catfile( $plot_dir, $merged_base . '_PCA.png');
		my $distance = catfile( $plot_dir, $merged_base . '_distance.png');
		my $pearson  = catfile( $plot_dir, $merged_base . '_pearson.png');
		$string .= <<END;

### Merged-replicate Correlation Metrics

Depth-normalized, scaled fragment counts were collected over each interval for direct
comparison analysis. The Principle Component Analysis (PCA) plot of the top two
components is plotted below, along with the Euclidean distance and Pearson
correlation plots.

![rep-merge_pca]($pca)

![rep-merge_distance]($distance)

![rep-merge_pearson]($pearson)
END
	}

	# Merged replicates profile
	if ( $self->plot and $merged_count ) {
		my $size     = ( 25 * $self->binsize ) / 1000;
		my $tot_size = $size * 2;
		my ( $hm_plot, $which_plot );
		if ($rep_num > 12) {
			$hm_plot   = catfile( $plot_dir, $merged_base .
				'_profile_mean_fragment_sorted_hm.png' );
			$which_plot = 'individual';
		}
		elsif ($rep_num >= 2 and $rep_num <= 12) {
			$hm_plot   = catfile( $plot_dir, $merged_base .
				'_profile_replicate_fragment_sorted_hm.png' );
			$which_plot = 'mean';
		}
		my $line_plot = catfile( $plot_dir, $merged_base . 
			'_profile_mean_fragment_summary.png' );
		if ( -e catfile( $self->dir, $hm_plot ) ) {
			$string .= <<END;

### Merged-replicate Fragment Profile

Fragment coverage data was collected over a $tot_size Kb window centered on each
interval midpoint (+/- $size Kb). Data was collected for both sample mean coverage
tracks and individual replicate coverage tracks. The $which_plot fragment coverage
heat map is shown below. The peaks (rows) are grouped based on which sample generated
the peak call and indicated on the left. The mean fragment class-average summary line
plot for all peaks is also shown.  Additional class-average summary line plots for
each group, as well as additional heat maps, are available in the **$plot_dir**
folder.

![rep-merge_replicate_frag_profile]($hm_plot)

![rep-merge_mean_frag_summary]($line_plot)
END
		}
	}

	# Enrichment plot
	if ( $self->plot and $merged_count ) {
		$string .= <<END;

### Merged-replicate Enrichment

The mean log2 Fold Enrichment value from the replicate-mean coverage was collected over
each peak interval.
END
		my $hm_plot    = catfile( $plot_dir, $merged_base . '_log2FE_hm.png' );
		my $hm_k4_plot = catfile( $plot_dir, $merged_base . '_log2FE_hm_K4.png' );
		if ( -e catfile( $self->dir, $hm_k4_plot ) ) {
			$string .= <<END;
The values were plotted as a heat map with a k-means clustering of 4. Additional 
higher k value plots are also available in the **$plot_dir** folder.

![rep-merge_enrichment]($hm_k4_plot)
END
		}
		elsif ( -e catfile( $self->dir, $hm_plot ) ) {
			$string .= <<END;
The values were plotted as a sorted heat map. Not enough intervals were available for 
k-means clustering.

![rep-merge_enrichment]($hm_plot)
END
		}
	}

	# Efficiency plot
	my $eff_plot = catfile( $plot_dir, $merged_base . '.chip_efficiency.png' );
	if ( $self->plot and -e catfile( $self->dir, $eff_plot) ) {
		$string .= <<END;

### Merged-replicate Efficiency

The Fraction of Reads (fragments) In Peaks, or FRIP, score is a measurement of
efficiency or specificity. The higher the fraction, the more specific or efficient
the antibody, targeting, or enrichment. It should be high for ChIP samples and low
for reference control (Input) samples. Fragment counts were collected over intervals
and expressed as a fraction of total fragments observed. The efficency plot is
plotted with countsfrom all replicates.

![rep-merge_efficiency]($eff_plot)
END
	}

	return $string;
}

sub add_merged_replicates_broad_report {
	my $self = shift;

	my $overlap = $self->minpeakover || 'n-1';
	my $gap     = $self->fragsize;
	my $rep_num = 0;  # total number of replicates
	my $string  = <<END;

--------

## Gapped Merged-replicate Peaks

Independently called replicate gapped peaks were merged for each sample, requiring at
least $overlap overlaps and a maximum gap of $gap bp. The final number of peaks for
each sample are below.

| Sample | File | Number |
|---|---|---|
END
	foreach my $Job ($self->list_jobs) {
		next unless $Job->repmerge_gappeak;
		$string .= sprintf "| %s | `%s` | %s |\n", $Job->job_name,
			( splitpath( $Job->repmerge_gappeak ) )[2],
			format_with_commas( $self->count_file_lines( $Job->repmerge_gappeak ) );
		$rep_num += scalar( $Job->chip_rep_names );
	}
	
	# Merged replicates
	my $merged_base  = ( splitpath( $self->repmerge_merge_base ) )[2];
	$merged_base    .= '_broad';
	my $merged_peak  = $merged_base . '.bed';
	my $plot_dir     = 'Replicate-Merge_Plots';
	my $merged_count = format_with_commas( $self->count_file_lines(
		$self->repmerge_merge_base . '_broad.bed' ) );
	$string .= <<END;

These merged-sample peaks were then merged into a final master broad replicate-merge
peak set with **$merged_count** intervals. This file, `$merged_peak`, was then used
for sample comparison.
END

	# add plots
	if ( $self->plot and $merged_count ) {
		my $merge_upset  = catfile( $plot_dir, $merged_base . '.intersection_upset.png' );
		my $merge_length = catfile( $plot_dir, $merged_base . '.peak_lengths.png' );
		$string .= <<END;
The UpSet intersection between samples and sample peak length distributions are shown
below. Additional plots may be found in the **$plot_dir** folder.

![gap_rep-merge_upset]($merge_upset)

![gap_rep-merge_peak_length]($merge_length)
END
	}

	# Merged replicates correlation
	if ( $self->plot and length($merged_count) > 1 and $rep_num > 2 ) {
		my $pca      = catfile( $plot_dir, $merged_base . '_PCA.png');
		my $distance = catfile( $plot_dir, $merged_base . '_distance.png');
		my $pearson  = catfile( $plot_dir, $merged_base . '_pearson.png');
		$string .= <<END;

### Gapped Merged-replicate Correlation Metrics

Depth-normalized, scaled fragment counts were collected over each broad interval for
direct comparison analysis. The Principle Component Analysis (PCA) plot of the top
two components is plotted below, along with the Euclidean distance and Pearson
correlation plots.

![gap_rep-merge_pca]($pca)

![gap_rep-merge_distance]($distance)

![gap_rep-merge_pearson]($pearson)
END
	}

	# Enrichment plot
	if ( $self->plot and $merged_count ) {
		$string .= <<END;

### Gapped Merged-replicate Enrichment

The mean log2 Fold Enrichment value from the replicate-mean coverage was collected over
each broad peak interval.
END
		my $hm_plot    = catfile( $plot_dir, $merged_base . '_log2FE_hm.png' );
		my $hm_k4_plot = catfile( $plot_dir, $merged_base . '_log2FE_hm_K4.png' );
		if ( -e catfile( $self->dir, $hm_k4_plot ) ) {
			$string .= <<END;
The values were plotted as a heat map with a k-means clustering of 4. Additional 
higher k value plots are also available in **$plot_dir**.

![gap_rep-merge_]($hm_k4_plot)
END
		}
		elsif ( -e catfile( $self->dir, $hm_plot ) ) {
			$string .= <<END;
The values were plotted as a sorted heat map. Not enough intervals were available for 
k-means clustering.

![gap_rep-merge_enrichment]($hm_plot)
END
		}
	}

	$string .= <<END;

Profile fragment or log2 fold enrichment plots are not generated for broad, gapped 
peaks due to the variably large size of the broad intervals.
END
	return $string;
}

sub add_mean_replicates_report {
	my $self = shift;

	# Mean replicates
	my $rep_num  = 0;  # total number of replicates
	my $peaksize = $self->peaksize;
	my $peakgap  = $self->peakgap;
	my $cutoff   = $self->cutoff;
	my $decimal;
	{
		my $number  = int($cutoff) + 1;
		my $pattern = "%.$number" . 'f';
		$decimal    = sprintf $pattern, ( 1 / ( 10 ** $cutoff ) );
	}
	my $string = <<END;

--------

## Replicate-mean Peaks

Peaks were called for each sample from the sample fragment coverage and enrichment
tracks as described above. If replicates were provided, they were averaged
together. The number of peaks identified are listed below.

| Sample | File | Number |
|---|---|---|
END

	foreach my $Job ( $self->list_jobs ) {
		next unless $Job->repmean_peak;
		$string .= sprintf "| %s | `%s` | %s |\n", $Job->job_name,
			( splitpath( $Job->repmean_peak ) )[2],
			format_with_commas( $self->count_file_lines ( $Job->repmean_peak ) );
		$rep_num += scalar( $Job->chip_rep_names );
	}

	my $mean_base  = $self->out;
	my $plot_dir   = 'Plots' . $self->dir_suffix;
	if ( $self->independent ) {
		$mean_base  = ( splitpath( $self->repmean_merge_base ) )[2];
		$plot_dir   = 'Replicate-Mean_Plots';
	}
	my $mean_peak  = $mean_base . '.bed';
	my $mean_count = format_with_commas( $self->count_file_lines( 
		catfile( $self->dir, $mean_peak ) ) );
	$string .= <<END;

These peaks were then merged into a final master replicate-mean peak set with
**$mean_count** intervals. This file, `$mean_peak`, was then used for sample
comparison.
END

	# add plots
	if ( $self->plot and $mean_count ) {
		my $mean_upset  = catfile( $plot_dir, $mean_base . '.intersection_upset.png' );
		my $mean_length = catfile( $plot_dir, $mean_base . '.peak_lengths.png' );
		$string .= <<END;
The UpSet intersection between samples and sample peak length distributions are shown
below. Additional plots may be found in the **$plot_dir** folder.

![rep-mean_upset]($mean_upset)

![rep-mean_peak_length]($mean_length)
END
	}

	# Mean replicates correlation
	if ( $self->plot and length($mean_count) > 1 and $rep_num > 2 ) {
		my $pca      = catfile( $plot_dir, $mean_base . '_PCA.png' );
		my $distance = catfile( $plot_dir, $mean_base . '_distance.png' );
		my $pearson  = catfile( $plot_dir, $mean_base . '_pearson.png' );
		$string .= <<END;

### Replicate-mean Correlation Metrics

Depth-normalized, scaled fragment counts were collected over each interval for direct
comparison analysis. The Principle Component Analysis (PCA) plot of the top two
components is plotted below, along with the Euclidean distance and Pearson
correlation plots.

![rep-mean_pca]($pca)

![rep-mean_distance]($distance)

![rep-mean_pearson]($pearson)
END
	}

	# Mean replicates profile
	if ( $self->plot and $mean_count ) {
		my $size     = ( 25 * $self->binsize ) / 1000;
		my $tot_size = $size * 2;
		my $hm_plot   = catfile( $plot_dir, $mean_base .
			'_profile_mean_fragment_sorted_hm.png' );
		my $line_plot = catfile( $plot_dir, $mean_base .
			'_profile_mean_fragment_summary.png' );
		$string .= <<END;

### Replicate-mean Fragment Profile

Fragment coverage data was collected over a $tot_size Kb window centered on each peak
interval midpoint (+/- $size Kb). The replicate-mean fragment coverage heat map is
shown below. The peaks (rows) are grouped based on which sample generated the peak
call and indicated on the left. The mean fragment class-average summary line plot for
all peaks is also shown. Additional class-average summary line plots for each group,
as well as additional heat maps, are available in the **$plot_dir** folder.

![rep-mean_replicate_frag_profile]($hm_plot)

![rep-mean_mean_frag_summary]($line_plot)
END
	}

	# Enrichment plot
	if ( $self->plot and $mean_count ) {
		$string .= <<END;

### Replicate-mean Mean Enrichment

The mean log2 Fold Enrichment value was collected over each peak interval.
END

		my $hm_plot    = catfile( $plot_dir, $mean_base . '_log2FE_hm.png' );
		my $hm_k4_plot = catfile( $plot_dir, $mean_base . '_log2FE_hm_K4.png' );
		if ( -e catfile( $self->dir, $hm_k4_plot ) ) {
			$string .= <<END;
The values were plotted as a heat map with a k-means clustering of 4. Additional
higher k value plots are also available in the **$plot_dir** folder.

![rep-mean_enrichment]($hm_k4_plot)
END
		}
		elsif ( -e catfile( $self->dir, $hm_plot ) ) {
			$string .= <<END;
The values were plotted as a sorted heat map. Not enough intervals were available for 
k-means clustering.

![rep-mean_enrichment]($hm_plot)
END
		}
		else {
			$string .= <<END;
Unable to plot this map.

END
		}
	}

	# Efficiency plot
	my $eff_plot = catfile( $plot_dir, $mean_base . '.chip_efficiency.png' );
	if ( $self->plot and -e catfile( $self->dir, $eff_plot ) ) {
		$string .= <<END;

### Replicate-mean Efficiency

The Fraction of Reads (fragments) In Peaks, or FRIP, score is a measurement of
efficiency or specificity. The higher the fraction, the more specific or efficient
the antibody, targeting, or enrichment. It should be high for ChIP samples and low
for reference control (Input) samples. Fragment counts were collected over intervals
and expressed as a fraction of total fragments observed. The efficency plot is
plotted with counts from all replicates.

![rep-mean_efficiency]($eff_plot)
END
	}

	return $string;
}

sub add_mean_replicates_broad_report {
	my $self = shift;

	my $string = <<END;

--------

## Gapped Replicate-mean Peaks

Broad, gapped peaks were called for each sample from the replicate-mean fragment
coverage and enrichment tracks. Significant peaks were called as described above.
The number of peaks identified are listed below.

| Sample | File | Number |
|---|---|---|
END

	# add gapped peak files
	my $rep_num  = 0;  # total number of replicates
	foreach my $Job ( $self->list_jobs ) {
		next unless $Job->repmean_gappeak;
		$string .= sprintf "| %s | `%s` | %s |\n", $Job->job_name,
			( splitpath( $Job->repmean_gappeak ) )[2],
			format_with_commas( $self->count_file_lines ( $Job->repmean_gappeak ) );
		$rep_num += scalar( $Job->chip_rep_names );
	}

	my $mean_base  = $self->out . '_broad';
	my $plot_dir   = 'Plots' . $self->dir_suffix;
	if ( $self->independent ) {
		$mean_base  = ( splitpath( $self->repmean_merge_base ) )[2];
		$mean_base .= '_broad';
		$plot_dir   = 'Replicate-Mean_Plots';
	}
	my $mean_peak  = $mean_base . '.bed';
	my $mean_count = format_with_commas( $self->count_file_lines( 
		catfile( $self->dir, $mean_peak ) ) );
	$string .= <<END;

These peaks were then merged into a final master broad replicate-mean peak set with
**$mean_count** intervals. This file, `$mean_peak`, was then used for sample
comparison.
END

	# add plots
	if ( $self->plot and $mean_count ) {
		my $mean_upset  = catfile( $plot_dir, $mean_base . '.intersection_upset.png' );
		my $mean_length = catfile( $plot_dir, $mean_base . '.peak_lengths.png' );
		$string .= <<END;
The UpSet intersection between samples and sample peak length distributions are shown
below. Additional plots may be found in the **$plot_dir** folder.

![gap_rep-mean_upset]($mean_upset)

![gap_rep-mean_peak_length]($mean_length)
END
	}

	# Mean replicates correlation
	if ( $self->plot and length($mean_count) > 1 and $rep_num > 2 ) {
		my $pca      = catfile( $plot_dir, $mean_base . '_PCA.png' );
		my $distance = catfile( $plot_dir, $mean_base . '_distance.png' );
		my $pearson  = catfile( $plot_dir, $mean_base . '_pearson.png' );
		$string .= <<END;

### Gapped Replicate-mean Correlation Metrics

Depth-normalized, scaled fragment counts were collected over each broad interval for
direct comparison analysis. The Principle Component Analysis (PCA) plot of the top
two components is plotted below, along with the Euclidean distance and Pearson
correlation plots.

![gap_rep-mean_pca]($pca)

![gap_rep-mean_distance]($distance)

![gap_rep-mean_pearson]($pearson)
END
	}

	# Enrichment plot
	if ( $self->plot and $mean_count ) {
		$string .= <<END;

### Gapped Replicate-mean Mean Enrichment

The mean log2 Fold Enrichment value was collected over each broad peak interval.
END

		my $hm_plot    = catfile( $plot_dir, $mean_base . '_log2FE_hm.png' );
		my $hm_k4_plot = catfile( $plot_dir, $mean_base . '_log2FE_hm_K4.png' );
		if ( -e catfile( $self->dir, $hm_k4_plot ) ) {
			$string .= <<END;
The values were plotted as a heat map with a k-means clustering of 4. Additional
higher k value plots are also available in the **$plot_dir**.

![gap_rep-mean_enrichment]($hm_k4_plot)
END
		}
		elsif ( -e catfile( $self->dir, $hm_plot ) ) {
			$string .= <<END;
The values were plotted as a sorted heat map. Not enough intervals were available for 
k-means clustering.

![gap_rep-mean_enrichment]($hm_plot)
END
		}
	}


	$string .= <<END;

Profile fragment or log2 fold enrichment plots are not generated for broad, gapped 
peaks due to the variably large size of the broad intervals.
END
	return $string;
}

sub add_mean_merge_comparison {
	my $self = shift;
	
	# get counts
	my $merged_peak  = ( splitpath( $self->repmerge_merge_base ) )[2] . '.bed';
	my $merged_count = format_with_commas( $self->count_file_lines(
		$self->repmerge_merge_base . '.bed' ) );
	my $mean_peak  = ( splitpath( $self->repmean_merge_base ) )[2] . '.bed';
	my $mean_count = format_with_commas( $self->count_file_lines( 
		catfile( $self->dir, $mean_peak ) ) );
	
	# get overlap stats
	my $intersect_file = catfile( $self->dir, 'Analysis' . $self->dir_suffix,
		'mean_merge.intersection.txt');
	my $Data = Bio::ToolBox->load_file($intersect_file) or
		do {
			carp "cannot read $intersect_file! $OS_ERROR";
			return;
		};
	my $jaccard       = $Data->value(2,3);
	my $count_number  = format_with_commas( $Data->value(2,4) );
	my $count_percent = sprintf "%.0f%%", $Data->value(2,5) * 100;
	my $evaluation;
	if ($jaccard > 0.9) {
		$evaluation = 'excellent';
	}
	elsif ($jaccard > 0.7) {
		$evaluation = 'good';
	}
	elsif ($jaccard > 0.5) {
		$evaluation = 'moderate';
	}
	else {
		$evaluation = 'poor';
	}
	$jaccard = sprintf "%.0f%%", $jaccard * 100;

	my $string = <<END;

--------

## Comparison between Replicate-Mean and Replicate-Merge Peaks

The two methods of calling peaks from replicates, Replicate-Mean and
Replicate-Merge, are comparable but mutually exclusive methods.
_Only one set should be chosen for further analysis._ They should not
be merged.

In ideal experiments, replicates that are strongly correlated should produce
nearly identical call sets between the two methods. Less correlative replicates,
or replicates with dissimilar strengths or efficiencies, may result in one
method working better than the other. 

Typically, if the replicates have good correlation, the Replicate-Merge peaks
should be used preferentially. Otherwise, the method generating the most peak
intervals should be used. 
END

	if ( $merged_count > $mean_count ) {
		$string .= <<END;
In this case, the **Replicate-Merge** method generated more peak intervals, with
$merged_count versus $mean_count intervals.
END
	}
	elsif ( $merged_count < $mean_count ) {
		$string .= <<END;
In this case, the **Replicate-Mean** method generated more peak intervals, with
$mean_count versus $merged_count intervals.
END
	}
	elsif ( $merged_count == $mean_count ) {
		$string .= <<END;
In this case, the two methods generated an identical number of peak intervals,
$mean_count, so the methods are equivalent.
END
	}

	$string .= <<END;

The spatial overlap (Jaccard) between the two methods is $jaccard, which is 
$evaluation. There were $count_number overlapping intervals ($count_percent).
END
	
	# add plots when available
	if ($self->plot) {
		my $plot_dir = 'Replicate-Merge_Plots';
		my $upset_image = catfile('Replicate-Merge_Plots' . $self->dir_suffix,
			'mean_merge.intersection_upset.png');
		my $spatial_image = catfile('Replicate-Merge_Plots' . $self->dir_suffix,
			'mean_merge.spatial_pie.png');
		$string .= <<END;

Below are the count intersection UpSet and spatial overlap pie plots for the
intersection between Replicate-Merge and Replicate-Mean peak intervals.

![mean-merge-upset]($upset_image)

![mean-merge-pie]($spatial_image)

END
	}

	if ($self->broad) {
		$string .= <<END;

No comparisons are made for the gapped peak intervals.

END
	}
	return $string;
}


sub add_summary_report {
	my $self = shift;
	
	my $string = <<END;

--------

## Summary Tables

Summary tables were written for the master list of peaks identified across all samples.
These include the coordinates, peak name, a boolean value (Y or N) indicating whether
the peak was identified in each sample, and the mean log2 Fold Enrichment and
maximum observered Q-Value statistical values over each interval for each sample.
These files are tab-delimited text files and can be directly imported into spreadsheet
programs.
END
	my ( $mean_count, $merge_count, $samples );
	if ( $self->independent ) {
		my $base       = ( splitpath( $self->repmerge_merge_base ) )[2];
		my $merge_file = $base . '_summary.tsv';
		$merge_count   = $base . '_counts.txt.gz';
		$base          = ( splitpath( $self->repmean_merge_base ) )[2];
		my $mean_file  = $base . '_summary.tsv';
		$mean_count    = $base . '_counts.txt.gz';
		$samples       = $base . '_samples.txt';
		$string .= <<END;
The merged-replicate peaks are in `$merge_file`.
The replicate-mean peaks are in `$mean_file`.
END
	}
	else {
		my $mean_file = $self->out . '_summary.tsv';
		$mean_count   = $self->out . '_counts.txt.gz';
		$samples      = $self->out . '_samples.txt';
		$string .= <<END;
The replicate-mean peaks are in `$mean_file`.
END
	}

	my $analdir = 'Analysis' . $self->dir_suffix;
	$string .= <<END;

Count tables are in the **$analdir** folder and may be used directly in differential
analysis between conditions, such as DESeq2. The counts are depth-normalized
and scaled and should be used directly without additional normalization (otherwise
any differences may be normalized away). For example, see the script `run_DESeq2.R`.
END
	if ( $self->independent ) {
		$string .= <<END;
The merged-replicate count table is `$merge_count`.
END
	}
	$string .= <<END;
The replicate-mean count table is `$mean_count`. The sample list table is `$samples`.

END
}

sub pandoc_header {
	my $self = shift;
	return <<END;
<style>
html {
  line-height: 1.5;
  font-family: Times, serif;
  font-size: 18px;
  color: #1a1a1a;
  background-color: #fdfdfd;
}
body {
  margin: 0 auto;
  max-width: 50em;
  padding-left: 50px;
  padding-right: 50px;
  padding-top: 50px;
  padding-bottom: 50px;
  hyphens: none;
  word-wrap: break-word;
  text-rendering: optimizeLegibility;
  font-kerning: normal;
}
code {
  font-family: Menlo, Monaco, Consolas, 'Lucida Console', monospace;
  font-size: 75%;
  margin: 0;
  hyphens: manual;
  word-wrap: break-word;
}
img {
  max-width: 75%;
}
table {
  margin: 1em 0;
  border-collapse: collapse;
  width: 100%;
  overflow-x: auto;
  display: block;
  font-variant-numeric: lining-nums tabular-nums;
}
tbody {
  margin-top: 0.5em;
  border-top: 1px solid #1a1a1a;
  border-bottom: 2px solid #1a1a1a;
}
th {
  border-top: 2px solid #1a1a1a;
  padding: 0.25em 0.5em 0.25em 0.5em;
}
td {
  padding: 0.125em 0.5em 0.25em 0.5em;
}
</style>
END
}

1;

=head1 name

Bio::MultiRepChIPSeq::report - Methods for generating markdown report

=head1 DESCRIPTION

Additional methods for the L<Bio::MultiRepChIPSeq::Runner> object for generating
a Markdown formatted report of the pipeline results. This uses GitHub-Flavored 
Markdown, primarily for table support.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

=head1 LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. 



