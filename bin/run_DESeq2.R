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
              help="Name of first ChIP condition"),
  make_option(c("-s","--second"), default="NA",
              help="Name of second ChIP condition or Input reference"),
  make_option(c("-t","--threshold"), default=0.01,
              help="Threshold adjusted p-value for filtering, default 0.01"),
  make_option(c("-m","--min"), default=50,
              help="Minimum base count sum, default 50"),
  make_option(c("--all"), action="store_true", default=FALSE,
              help="Report all windows, not just significant")
)

parser <- OptionParser(option_list=opts, description = "
  This script will run a basic DESeq2 differential analysis 
  between two conditions to identify differential (or enriched 
  in the case of ChIP and Input) ChIPseq regions. It requires 
  an input text file with chromosome, start,and stop columns 
  (or a coordinate name string), along with columns of alignment 
  counts for both ChIP1 and ChIP2 (or Input reference) sample 
  replicates. DESeq2 requires ideally 3 or more replicates per 
  condition to estimate variance for normalization and significance 
  testing. A sample condition file is required, consisting of two 
  columns, sample identifiers and groups (conditions).
  
  Result intervals are filtered for the given adjusted  
  threshold as well as a minimum base count (mean of ChIP1 and 
  ChIP2 replicate counts). Results are written with A bedGraph file of regularized, log2 
  differences between ChIP1 and ChIP2 (or Reference) is written 
  for converting into a bigWig for visualization. Merged 
  significant intervals of enrichment are written as a bed file.
  
")

opt <- parse_args(parser)
if( opt$count == "NA" ){
  print_help(parser)
  quit(status=1)
}

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(DESeq2))


# input sample file
samples <- read.delim(opt$sample, header = TRUE, row.names=1, 
                      check.names = F, sep = "\t", comment.char = "#")
idx <- which(samples[,1] == opt$first | samples[,1] == opt$second)
cond <- data.frame(row.names = rownames(samples)[idx], condition = samples[idx,1])

# input count file
counts <- read.delim(opt$count, header = TRUE, sep = "\t", 
                     check.names = F, comment.char = "#", na.strings = ".")
counts[is.na(counts)] <- 0

# output basename
if (opt$output == "NA"){
  opt$output <- paste0(opt$first,'_',opt$second)
}


# build genomic ranges for later
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


# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData=counts[,rownames(cond)], 
                              colData=cond, 
                              design =~ condition)
dds <- DESeq(dds)
ddsResults1 <- results(dds, alpha = opt$threshold, contrast = c('condition', opt$first, opt$second))

# shrink log fold changes - is this necessary?
# ddsResults1 <- lfcShrink(dds=dds, res=ddsResults1, contrast= c('condition', opt$first, opt$second))


# write all results
if (opt$all == TRUE) {
	allresults <- data.frame(
	  Chromosome = seqnames(targets),
	  Start = start(targets),
	  End = end(targets),
	  baseMean = ddsResults1$baseMean,
	  Log2FoldChange = ddsResults1$log2FoldChange,
	  Padj = ddsResults1$padj
	)
    # check for name and add it if present
    if (length(namecol)) {
      allresults$Name <- counts[,namecol]
      allresults <- allresults[,c(1,2,3,7,4:6)]
    }
    # write results
    write.table(allresults,file=paste0(opt$output,".all.txt"),
                sep = "\t", row.names = F, quote = F, col.names = T)
}




# write bedGraph of log fold changes
write.table(data.frame(seqnames(targets),start(targets),end(targets), 
                       round(ddsResults1$log2FoldChange,3)),
            file=paste0(opt$output,"_logEnrichment.bdg"),
            sep = "\t", row.names = F, quote = F, col.names = F)


# select significant results
idx <- which(ddsResults1$padj <= opt$threshold & ddsResults1$baseMean >= opt$min)
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

# write first ChIP enriched merged bed file
idx <- which(results$Log2FoldChange > 0)
sig2.gr <- reduce(sig.gr[idx])
if (length(sig2.gr)) {
  write.table(data.frame(seqnames(sig2.gr),start(sig2.gr),end(sig2.gr)), 
              file=paste0(opt$output, '_', opt$first, 'Only.bed'),
              sep = "\t", row.names = F, quote = F, col.names = F)
}


# write second ChIP enriched merged bed file
idx <- which(results$Log2FoldChange < 0)
sig3.gr <- reduce(sig.gr[idx])
if (length(sig3.gr)) {
  write.table(data.frame(seqnames(sig3.gr),start(sig3.gr),end(sig3.gr)), 
              file=paste0(opt$output, '_', opt$second, 'Only.bed'),
              sep = "\t", row.names = F, quote = F, col.names = F)
}


