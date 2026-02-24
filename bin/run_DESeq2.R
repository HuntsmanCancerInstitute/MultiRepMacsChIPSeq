#!/usr/bin/env Rscript

# script for running a basic DESeq2 differential ChIPseq analysis

# Timothy J. Parnell, PhD
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112
#  
# This package is free software; you can redistribute it and/or modify
# it under the terms of the Artistic License 2.0.  
# 
# Updated versions of this file may be found in the repository
# https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq


suppressPackageStartupMessages(library("optparse"))

opts <-  list(
  make_option(c("-c", "--count"), default="NA",
              help="Input file containing count data"),
  make_option(c("-a", "--sample"), default="NA",
              help="Sample condition file containing identifiers and conditions"),
  make_option(c("-o", "--output"), default="NA",
              help="Output file basename, default 'first_second' names"),
  make_option(c("-f","--first"), default="NA",
              help="Name of first ChIP condition, numerator"),
  make_option(c("-s","--second"), default="NA",
              help="Name of second ChIP condition or Input reference, denominator"),
  make_option(c("-b","--batch"), action="store_true", default=FALSE,
              help="Use third column in customized sample table as batch factor"),
  make_option(c("-t","--threshold"), default=0.1,
              help="Threshold adjusted p-value (alpha) for significance (FDR), default 0.1"),
  make_option("--norm", action="store_true", default=FALSE,
              help="Input counts are already normalized"),
  make_option("--all", action="store_true", default=FALSE,
              help="Report all windows, not just significant"),
  make_option("--shrink", action="store_true", default=FALSE,
              help="Shrink Log2FC values for low count peaks")
)

parser <- OptionParser(option_list=opts, description = "
  This script will run a basic DESeq2 differential analysis 
  between two conditions to identify differential (or enriched 
  in the case of ChIP and Input) ChIPseq regions. 
  
  It requires an input text file with chromosome, start, and stop 
  columns (or a coordinate name string), along with columns of fragment 
  counts for both ChIP1 and ChIP2 (or Input reference) sample replicates.
  DESeq2 requires at least 2 (ideally more) replicates per condition 
  to estimate variance for normalization and significance testing.
  
  Count data may be already normalized to a uniform genomic depth, in
  which case it should indicated as such and SizeFactors will be set to 1.
  Otherwise counts will be normalized by DESeq2. WARNING: Default 
  normalization by DESeq2 may potentially erase any enrichment signal.
  
  A sample condition file is required, consisting of two columns, 
  unique sample identifiers and group names (conditions). If desired,
  a third column could be added as an additional batch factor. Only
  two groups are used in the contrast, defined with the --first and 
  --second options; any additional groups and samples are ignored.
  
  Results are filtered at the indicated threshold or alpha level 
  (False Discovery Rate). Significant regions are written to a text file
  and respective bed files. Log2 Fold Changes may be shrunken to reduce the
  effect of low count peaks. If additional filtering on Log2FC will be 
  performed, the values should be shrunk.
  
")

opt <- parse_args(parser)
if( opt$count == "NA" ){
  print_help(parser)
  quit(status=1)
}

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(DESeq2))


# input sample file
# only first 2 or 3 columns are used
samples <- read.delim(opt$sample, header = TRUE, row.names=1, 
                      check.names = F, sep = "\t", comment.char = "#")
idx <- which(samples[,1] == opt$first | samples[,1] == opt$second)
cond <- data.frame(row.names = rownames(samples)[idx], condition = samples[idx,1])

# use 3rd column as batch control
if (opt$batch == TRUE) {
  cond$batch <- samples[idx,2]
}

# input count file
counts <- read.delim(opt$count, header = TRUE, sep = "\t", 
                     check.names = F, comment.char = "#", na.strings = ".")
counts[is.na(counts)] <- 0

# output basename
if (opt$output == "NA"){
  opt$output <- paste0(opt$first,'_',opt$second)
}


# load count table
# build genomic ranges for convenient use later - though probably not necessary
chromcol <- grep("^Chromo",colnames(counts))
startcol <- grep("^Start",colnames(counts))
stopcol  <- grep("^Stop|End", colnames(counts))
namecol  <- grep("^Name", colnames(counts))
idcol  <- grep("^Primary_ID", colnames(counts))
if (length(chromcol) & length(startcol) & length(stopcol)) {
  # this is easy
  targets <- GRanges(seqnames=counts[,chromcol],
                     ranges=IRanges(counts[,startcol],counts[,stopcol]),
                     strand=NULL)
} else if (length(idcol)) {
  # must extract from id, hope it's as chr1:123-456 format
  coord <- do.call(rbind, strsplit(as.character(counts[,idcol]), ":|-"))
  targets <- GRanges(seqnames=coord[,1],
                     ranges=IRanges(as.integer(coord[,2]),as.integer(coord[,3])),
                     strand=NULL)
} else if (length(namecol)) {
  # try extracting from name, hope it's as chr1:123-456 format
  coord <- do.call(rbind, strsplit(as.character(counts[,namecol]), ":|-"))
  targets <- GRanges(seqnames=coord[,1],
                     ranges=IRanges(as.integer(coord[,2]),as.integer(coord[,3])),
                     strand=NULL)
}
message("Loaded ", nrow(counts), " peaks for analysis between ", opt$first, " and ", opt$second)

# set design formula
if (opt$batch == TRUE) {
  des <- formula( ~ batch + condition)
} else {
  des <- formula( ~ condition)
}

# Run DESeq2 using only the indicated samples from the count table
dds <- DESeqDataSetFromMatrix(countData=counts[,rownames(cond)], 
                              colData=cond, 
                              design = des)

# default DESeq2 normalization can destroy biological significance
# since peak counts represent a tiny fraction of actual sequencing depth
# best to normalize to sequencing depth during counting step elsewhere
if (opt$norm == TRUE) {
  sizeFactors(dds) <- rep(1, nrow(cond))
}

# generate results using local fitting since it seems to work better
dds <- DESeq(dds, fitType = 'local')
ddsResults1 <- results(dds, alpha = opt$threshold, contrast = c('condition', opt$first, opt$second))
message(summary(ddsResults1))

# shrink log fold changes
# using the original normal method as the default apeglm is too aggressive here
if (opt$shrink == TRUE) {
  ddsResults1 <- lfcShrink(dds=dds, res=ddsResults1, type="normal", quiet=TRUE,
                           contrast= c('condition', opt$first, opt$second))
}


# write all results
if (opt$all == TRUE) {
	allresults <- data.frame(
	  Chromosome = seqnames(targets),
	  Start = start(targets),
	  End = end(targets),
	  baseMean = ddsResults1$baseMean,
	  Log2FoldChange = ddsResults1$log2FoldChange,
	  PValue = ddsResults1$pvalue,
	  Padj = ddsResults1$padj
	)
    # check for name and add it if present
    if (length(namecol)) {
      allresults$Name <- counts[,namecol]
      allresults <- allresults[,c(1,2,3,8,4:7)]
    }
    # write results
    write.table(allresults,file=paste0(opt$output,".all.txt"),
                sep = "\t", row.names = F, quote = F, col.names = T)
    message("wrote all results to file ", opt$output, ".all.txt")
}




# write bedGraph of log fold changes
# write.table(data.frame(seqnames(targets),start(targets),end(targets), 
#                        round(ddsResults1$log2FoldChange,3)),
#             file=paste0(opt$output,"_logEnrichment.bdg"),
#             sep = "\t", row.names = F, quote = F, col.names = F)


# select significant results
idx <- which(ddsResults1$padj <= opt$threshold)
sig.gr <- targets[idx]
results <- data.frame(
  Chromosome = seqnames(sig.gr),
  Start = start(sig.gr),
  End = end(sig.gr),
  baseMean = ddsResults1$baseMean[idx],
  Log2FoldChange = ddsResults1$log2FoldChange[idx],
  Padj = ddsResults1$padj[idx]
)

# check for name and add it back if present
if (length(namecol)) {
  results$Name <- counts[idx,namecol]
  results <- results[,c(1,2,3,7,4:6)]
}

# add normalized sample counts
results <- cbind(results, counts(dds,normalized=TRUE)[idx,])

# write results
write.table(results,file=paste0(opt$output,".txt"),
            sep = "\t", row.names = F, quote = F, col.names = T)
message("wrote ", nrow(results), " significant results to file ", opt$output, ".txt")


# write first ChIP enriched merged bed file
bed1 <- results[results$Log2FoldChange > 0,]
if (nrow(bed1)) {
  if (length(namecol)) {
    bed1 <- bed1[,c(1:4,7)]
  } else {
    bed1 <- bed1[,c(1:3,6)]
    bed1$Name <- rep(opt$first, times = nrow(bed1))
    bed1 <- bed1[,c(1:3,5,4)]
  }
  bed1$Padj <- round(-log10(bed1$Padj), digits = 2)
  write.table(bed1, file=paste0(opt$output, "_", opt$first, "Only.bed"),
              sep = "\t", row.names = F, quote = F, col.names = F)
  message("wrote ", nrow(bed1), " ", opt$first, " enriched regions to ", opt$output, "_",
          opt$first, "Only.bed")
}


# write second ChIP enriched merged bed file
bed2 <- results[results$Log2FoldChange < 0,]
if (nrow(bed2)) {
  if (length(namecol)) {
    bed2 <- bed2[,c(1:4,7)]
  } else {
    bed2 <- bed2[,c(1:3,6)]
    bed2$Name <- rep(opt$second, times = nrow(bed2))
    bed2 <- bed2[,c(1:3,5,4)]
  }
  bed2$Padj <- round(-log10(bed2$Padj), digits = 2)
  write.table(bed2, file=paste0(opt$output, "_", opt$second, "Only.bed"),
              sep = "\t", row.names = F, quote = F, col.names = F)
  message("wrote ", nrow(bed2), " ", opt$second, " enriched regions to ", opt$output, "_",
          opt$second, "Only.bed")
}

# plot volcano
library(ggplot2)
vdata <- data.frame(log2fc = ddsResults1$log2FoldChange, padj = ddsResults1$padj)
vdata <- na.omit(vdata)
fc <- max(abs(vdata$log2fc), na.rm=TRUE)
p1 <- ggplot(data = vdata, aes(x=log2fc, y=-log10(padj))) +
  geom_point(aes(fill = log2fc), color="gray20", shape = 21, size=1) +
  xlim(-fc, fc) + 
  theme_light() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-fc, fc), guide = "none") +
  ggtitle(paste0(opt$first, " vs ", opt$second))
ggsave(plot = p1, filename = paste0(opt$output, ".png"), 
       height = 6, width = 6)
message("wrote volcano plot to file ", opt$output, ".png")
