#!/usr/bin/env Rscript

# script for running a basic enrichR function from normR

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
  make_option(c("-i", "--input"), default="NA",
              help="Input file containing count data"),
  make_option(c("-o", "--output"), default="NA",
              help="Output file basename, default 'chip_ref' names"),
  make_option(c("-c","--chip"), default="NA",
              help="Name of ChIP count column"),
  make_option(c("-r","--ref"), default="NA",
              help="Name of reference (or Input) count column"),
  make_option(c("-t","--threshold"), default=0.001,
              help="Threshold Q-value for filtering, default 0.001"),
  make_option(c("-m","--min"), default=50,
              help="Minimum interval count sum, default 50"),
  make_option(c("--all"), action="store_true", default=FALSE,
              help="Report all windows, not just significant")
)

parser <- OptionParser(option_list=opts, description = "
 This script will run a basic enrichR function from the normR 
 package to identify enriched ChIPseq regions. It requires 
 an input text file with chromosome, start,and stop columns, 
 along with a columns of alignment counts for both ChIP and 
 reference control (Input) samples. 
 
 Result intervals are filtered for the given Q-value FDR 
 threshold, as well as a minimum read count (ChIP+Reference). 
 A bedGraph file of standardized, log enrichment between 
 ChIP and Reference is written for converting into a bigWig 
 for visualization. Merged significant intervals for 
 enrichment are written as bed files.
")

opt <- parse_args(parser)
if( opt$input == "NA" ){
  print_help(parser)
  quit(status=1)
}

suppressPackageStartupMessages(library(normr))
suppressPackageStartupMessages(library(GenomicRanges))

# input file
counts <- read.delim(opt$input, header = TRUE, 
                     check.names = F, sep = "\t", comment.char = "#")

# output basename
if (opt$output == "NA"){
  opt$output <- paste0(opt$chip,'_',opt$ref)
}

# columns
chromcol <- grep("^Chromo",colnames(counts))
startcol <- grep("^Start",colnames(counts))
stopcol  <- grep("^Stop|End", colnames(counts))
chipcol  <- grep(paste0('^',opt$chip,'$'), colnames(counts))
refcol   <- grep(paste0('^',opt$ref,'$'), colnames(counts))
namecol  <- grep("^Name", colnames(counts))


# target genomic ranges
if (length(chromcol) & length(startcol) & length(stopcol)) {
  # this is easy
  targets <- GRanges(seqnames=counts[,chromcol],
                     ranges=IRanges(counts[,startcol],counts[,stopcol]),
                     strand=NULL)
} else {
  # must extract from name, hope it's as chr1:123-456 format
  if (length(namecol)) {
    coord <- do.call(rbind, strsplit(as.character(counts[,namecol]), ":|-"))
    targets <- GRanges(seqnames=coord[,1],
                       ranges=IRanges(as.integer(coord[,2]),as.integer(coord[,3])),
                       strand=NULL)
  }
}

# check target widths
# the normR model fit assumes all target intervals are the same size
w <- width(targets[1])
if (all(width(targets)) != w) {
  warning("WARNING: tested target intervals are not identical in width! This may invalidate normR fit")
}


# run enrichment
e1 <- enrichR(treatment = counts[,chipcol],
              control = counts[,refcol],
              genome = targets,
              iterations = 10)

# write bedGraph
write.table(data.frame(seqnames(targets),start(targets),end(targets),getEnrichment(e1)),
            file=paste0(opt$output,"_logEnrichment.bdg"),
            sep = "\t", row.names = F, quote = F, col.names = F)

# write all results
if (opt$all == TRUE){
    allresults <- data.frame(
      Chromosome = seqnames(targets),
      Start = start(targets),
      End = end(targets),
      ChIPCount = getCounts(e1)$treatment,
      ControlCount = getCounts(e1)$control,
      Enrichment = getEnrichment(e1),
      Pvalue = getPvalues(e1),
      QValue = getQvalues(e1)
    )
    # check for name and add it if present
    if (length(namecol)) {
      allresults$Name <- counts[,namecol]
      allresults <- allresults[,c(1,2,3,9,4:8)]
    }
    # write results
    write.table(allresults,file=paste0(opt$output,".all.txt"),
                sep = "\t", row.names = F, quote = F, col.names = T)

}
 


# get significant results
idx <- which(getQvalues(e1) < opt$threshold & rowSums(do.call(cbind, getCounts(e1))) >= opt$min)
sig.gr <- getRanges(e1)[idx]
results <- data.frame(
  Chromosome = seqnames(sig.gr),
  Start = start(sig.gr),
  End = end(sig.gr),
  ChIPCount = getCounts(e1)$treatment[idx],
  ControlCount = getCounts(e1)$control[idx],
  Enrichment = getEnrichment(e1)[idx],
  Pvalue = getPvalues(e1)[idx],
  QValue = getQvalues(e1)[idx]
)

# check for name and add it if present
if (length(namecol)) {
  results$Name <- counts[idx,namecol]
  results <- results[,c(1,2,3,9,4:8)]
}

# write results
write.table(results,file=paste0(opt$output,".txt"),
            sep = "\t", row.names = F, quote = F, col.names = T)

# write merged bed file
sig2.gr <- reduce(sig.gr)
write.table(data.frame(seqnames(sig2.gr),start(sig2.gr),end(sig2.gr)), 
            file=paste0(opt$output,".bed"),
            sep = "\t", row.names = F, quote = F, col.names = F)


