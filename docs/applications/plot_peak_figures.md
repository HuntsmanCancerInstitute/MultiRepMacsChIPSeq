# MultiRepMacsChIPSeq - plot\_peak\_figures

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## plot\_peak\_figures.R

This script generates a number of heat maps and plots for identified peak calls, 
including the following:

- Stacked bar chart of unique and duplicate fragment counts
- Scatter plot of before and after duplication counts versus non-duplicate
- PCA plot between all sample replicates based on fragment counts in peak
- Multiple pairwise correlation (Pearson, Spearman, and Euclidean distance) 
  heat maps and clusters of the counts between sample replicates
- Bar chart plots for peak number and genomic space for each sample peaks
- Box and whisker plots for peak length distributions
- Line plot for ranked peak signals for each sample replicate
- Heat map and cluster of the number of peak intersections
- Heat map and cluster of jaccard (spatial overlap) statistic between peaks
- Pie chart of spatial overlap fraction of total merged peak coverage for 
  each sample
- UpSet plot of peak interactions
- UpSet plot of peak spatial overlaps (optional, compute intensive)
- Bar chart of the fraction of fragment counts in corresponding peaks 
  for each replicate, a measure of ChIP efficiency
- Heat map of the mean q-value scores for each ChIP over merged peaks
- Heat map of the mean log2 fold enrichment for each ChIP over merged peaks
  with k-means clustering (4, 6, and 8 clusters)
- Profile heat map of the fragment density over the midpoint of merged
  peaks for all samples, with and without k-means (4) clustering, or
  clustered by sample calls
- Profile heat map of the log2 Fold Enrichment over the midpoint of merged
  peaks for all samples, with and without k-means (4) clustering, or
  clustered by sample calls
- Mean profile line plot of fragment density over merged peaks, including
  for each group sample calls
- Mean profile line plot of log2 Fold Enrichment over merged peaks,
  including for each group of sample calls

Samples are identified by color palette: Try `Set1`, `Set2`, `Set3`, `Spectral`,
`Dark2`, or any other named palette in `RColorBrewer`. Note that excessive sample
numbers may exceed the number the colors in a given palette (8 to 12). 

USAGE: 

	plot_peak_figures.R [options]

OPTIONS:

	-i INPUT, --input=INPUT
		Path and basename to the multirep_macs2_pipeline combined output
	
	-n MIN, --min=MIN
		Minimum log2 Fold Enrichment value to plot, default -4
	
	-x MAX, --max=MAX
		Maximum log2 Fold Enrichment value to plot, default 4
	
	-q QMAX, --qmax=QMAX
		Maximum q-value to plot, default 30
	
	-r FMAX, --fmax=FMAX
		Maximum fragment value to plot, default 4
	
	--mincluster=MINCLUSTER
		Minimum number of peaks to generate k-means plots, default 500
	
	--spatial
	        Generate spatial UpSet plots (compute intensive), default F
	
	-p PALETTE --palette=PALETTE
		RColorBrewer palette for samples, default Set1
	
	-f FORMAT, --format=FORMAT
		Format of output file: png, pdf, default png
	
	-h, --help
		Show this help message and exit



