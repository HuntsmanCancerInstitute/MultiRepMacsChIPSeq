#!/usr/bin/env Rscript

# script for running a basic diffR function from normR

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
              help="Output file basename, default 'first_second' names"),
  make_option(c("-f","--first"), default="NA",
              help="Name of first ChIP count column"),
  make_option(c("-s","--second"), default="NA",
              help="Name of second ChIP count column"),
  make_option(c("-t","--threshold"), default=0.001,
              help="Threshold Q-value for filtering, default 0.001"),
  make_option(c("-m","--min"), default=50,
              help="Minimum interval count sum, default 50"),
  make_option(c("--all"), action="store_true", default=FALSE,
              help="Report all windows, not just significant")
)

parser <- OptionParser(option_list=opts, description = "
 This script will run a basic diffR function from the normR 
 package to identify differential ChIPseq regions. It requires 
 an input text file with chromosome, start,and stop columns, 
 along with a columns of alignment counts for both ChIP1 and 
 ChIP2 samples; input is not needed. 
 
 Result intervals are filtered for the given Q-value FDR 
 threshold as well as a minimum readcount (ChIP1+ChIP2). 
 A bedGraph file of standardized, log difference between 
 ChIP1 and ChIP2 is written for converting into a bigWig 
 for visualization. Merged significant intervals for 
 differential enrichment of ChIP1 (class 2) and ChIP2 
 (class 1) are written as bed files.
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
  opt$output <- paste0(opt$first,'_',opt$second)
}

# columns
chromcol <- grep("^Chromo",colnames(counts))
startcol <- grep("^Start",colnames(counts))
stopcol  <- grep("^Stop|End", colnames(counts))
chip1col  <- grep(paste0('^',opt$first,'$'), colnames(counts))
chip2col   <- grep(paste0('^',opt$second,'$'), colnames(counts))
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
d1 <- diffR(treatment = counts[,chip1col],
              control = counts[,chip2col],
              genome = targets,
              iterations = 10)

# write bedGraph
write.table(data.frame(seqnames(targets),start(targets),end(targets),
                       round(getEnrichment(d1),3)),
            file=paste0(opt$output,"_logDifference.bdg"),
            sep = "\t", row.names = F, quote = F, col.names = F)

# write all results
if (opt$all == TRUE){
    allresults <- data.frame(
      Chromosome = seqnames(targets),
      Start = start(targets),
      End = end(targets),
      ChIP1Count = getCounts(d1)$treatment,
      ChIP2Count = getCounts(d1)$control,
      Enrichment = getEnrichment(d1),
      Pvalue = getPvalues(d1),
      QValue = getQvalues(d1)
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
idx <- which(getQvalues(d1) < opt$threshold & 
               getClasses(d1) != "NA" & 
               rowSums(do.call(cbind, getCounts(d1))) >= opt$min)
sig.gr <- getRanges(d1)[idx]
results <- data.frame(
  Chromosome = seqnames(sig.gr),
  Start = start(sig.gr),
  End = end(sig.gr),
  Class = getClasses(d1)[idx],
  ChIP1Count = getCounts(d1)$treatment[idx],
  ChIP2Count = getCounts(d1)$control[idx],
  Enrichment = getEnrichment(d1)[idx],
  Pvalue = getPvalues(d1)[idx],
  QValue = getQvalues(d1)[idx]
)

# check for name and add it if present
if (length(namecol)) {
  results$Name <- counts[idx,namecol]
  results <- results[,c(1,2,3,10,4:9)]
}

# write results
write.table(results,file=paste0(opt$output,".txt"),
            sep = "\t", row.names = F, quote = F, col.names = T)


# write class 1 merged bed file
idx <- which(getQvalues(d1) <= opt$threshold & 
               getClasses(d1) == 1 & 
               rowSums(do.call(cbind, getCounts(d1))) >= opt$min)
if (length(idx)) {
  class1.gr <- reduce(getRanges(d1)[idx])
  write.table(data.frame(seqnames(class1.gr),start(class1.gr),end(class1.gr)), 
              file=paste0(opt$output,'_',opt$second,'Enriched.bed'),
              sep = "\t", row.names = F, quote = F, col.names = F)
}


# write class 2 merged bed file
idx <- which(getQvalues(d1) <= opt$threshold & 
               getClasses(d1) == 2 & 
               rowSums(do.call(cbind, getCounts(d1))) >= opt$min)
if (length(idx)) {
  class2.gr <- reduce(getRanges(d1)[idx])
  write.table(data.frame(seqnames(class2.gr),start(class2.gr),end(class2.gr)), 
              file=paste0(opt$output,'_',opt$first,'Enriched.bed'),
              sep = "\t", row.names = F, quote = F, col.names = F)
}



