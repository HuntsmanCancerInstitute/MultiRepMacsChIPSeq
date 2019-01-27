#!/usr/bin/env Rscript

# script for plotting collected peak information

suppressPackageStartupMessages(library("optparse"))

opts <-  list(
    make_option(c("-i", "--input"), default="NA",
         help="Path and basename to the multirep_macs2_pipeline combined output"),
    make_option(c("-n","--min"), default=-4,
         help="Minimum log2 Fold Enrichment value to plot, default -4"),
    make_option(c("-x","--max"), default=4,
         help="Maximum log2 Fold Enrichment value to plot, default 4"),
    make_option(c("-q","--qmax"), default=30,
         help="Maximum q-value to plot, default 30"),
    make_option(c("-f","--format"), default="png",
         help="Format of output file: png, pdf, default png")
)

parser <- OptionParser(option_list=opts, description = "
This script generates a number of heat maps and plots for identified peak calls, 
including the following:
  - heat map and cluster of jaccard (spatial overlap) statistic between peaks
  - heat map and cluster of the number of peak intersctions
  - heat map of the mean q-value scores for each ChIP over merged peaks
  - heat map of the mean log2 fold enrichment for each ChIP over merged peaks
    with k-means clustering (6, 8, and 10 clusters)
  - pairwise Pearson, Spearman, and Euclidean distance between all sample 
    replicates with heat map and cluster
  - PCA plot between all sample replicates
")

opt <- parse_args(parser)
if( opt$input == "NA" ){
 print_help(parser)
 quit(status=1)
}


library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# annotation colors for different numbers of clusters
# http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=4
color4 = c("1" = "#a6611a", "2" = "#dfc27d", "3" = "#80cdc1", "4" = "#018571")
color6 = c("1" = "#8c510a", "2" = "#d8b365", "3" = "#f6e8c3",
           "4" = "#c7eae5", "5" = "#5ab4ac", "6" = "#01665e")
color8 = c("1" = "#8c510a", "2" = "#bf812d", "3" = "#dfc27d",
           "4" = "#f6e8c3", "5" = "#c7eae5", "6" = "#80cdc1",
           "7" = "#35978f", "8" = "#01665e")
color10 = c("1" = "#543005", "2" = "#8c510a", "3" = "#bf812d",
            "4" = "#dfc27d", "5" = "#f6e8c3", "6" = "#c7eae5",
            "7" = "#80cdc1", "8" = "#35978f", "9" = "#01665e",
            "10" = "#003c30")


#### Functions

makekmeans <-function(rdata, hm_min, hm_max, k, hm_color, rowColor, figmain, outbase) {
  set.seed(100)
  n <- ncol(rdata)
  o <- apply(rdata, 1, mean)
  rdata <- rdata[order(o, decreasing = T),]
  rdata <- replace(rdata, rdata > hm_max, hm_max)
  rdata <- replace(rdata, rdata < hm_min, hm_min)
  kresult <- kmeans(rdata, k, nstart = 100, iter.max = 100)
  kdata <- cbind(rdata, kresult$cluster)
  o <- order(kdata[, n+1])
  kdata <- kdata[o, ]
  rowAnnot <- data.frame(row.names = rownames(kdata),
                         "Cluster"=as.character(kdata$`kresult$cluster`))
  pheatmap(kdata[,c(1:n)], cluster_cols = T, color = hm_color, 
           border_color = NA, breaks = seq(hm_min,hm_max,length=256), 
           cluster_rows = F, show_rownames = F, 
           show_colnames = T, annotation_row = rowAnnot, 
           annotation_colors = list("Cluster" = rowColor), main = figmain,
           filename = paste0(outbase, ".", opt$format), width = 8, height = 10)
  write.table(data.frame(kresult$cluster), file = paste0(outbase, "_klist.txt"), 
              sep = "\t", quote = F)
}

makehm <-function(rdata, hm_min, hm_max, hm_color, figmain, outbase) {
  o <- apply(rdata, 1, mean)
  rdata <- rdata[order(o, decreasing = T),]
  rdata <- replace(rdata, rdata > hm_max, hm_max)
  rdata <- replace(rdata, rdata < hm_min, hm_min)
  pheatmap(rdata, cluster_cols = T, color = hm_color, 
           border_color = NA, breaks = seq(hm_min,hm_max,length=256), 
           cluster_rows = F, show_rownames = F, show_colnames = T, main = figmain,
           filename = paste0(outbase, ".", opt$format), width = 8, height = 10)
}

makeDist <- function(rdata, cdata, outfile) {
  sampleDists = dist(t(rdata))
  sampleDistMatrix <- as.matrix( sampleDists )
  colours = colorRampPalette( rev(brewer.pal(9, 'Blues')) )(255)
  pheatmap(sampleDistMatrix, col=colours, annotation_row = cdata,
           cluster_cols = F, cluster_rows = T,
           main = "Distance Correlation of Sample Counts",
           filename = outfile, width = 10, height = 8)
}

makePearCorr <- function(rdata, cdata, outfile) {
  sampleCorr = cor(rdata, method = "pearson")
  sampleMatrix <- as.matrix( sampleCorr )
  colours = colorRampPalette( brewer.pal(9, 'Blues') )(255)
  pheatmap(sampleMatrix, col=colours, annotation_row = cdata, 
           cluster_cols = F, cluster_rows = T,
           main = "Pearson Correlation of Sample Counts",
           filename = outfile, width = 10, height = 8)
}

makeSpearCorr <- function(rdata, cdata, outfile) {
  sampleCorr = cor(rdata, method = "spearman")
  sampleMatrix <- as.matrix( sampleCorr )
  colours = colorRampPalette( brewer.pal(9, 'Blues') )(255)
  pheatmap(sampleMatrix, col=colours, annotation_row = cdata,
           cluster_cols = F, cluster_rows = T,
           main = "Spearman Correlation of Sample Counts",
           filename = outfile, width = 10, height = 8)
}


makePCA <- function(rdata, cdata, outfile) {
  if (ncol(rdata) != nrow(cdata)) {
    cdata <- subset(cdata, row.names(cdata) %in% colnames(rdata))
  }
  n  <- apply(rdata, 1, stats::var)
  x <- head(rdata[order(n, decreasing = T),], 1000)
  p <- prcomp(t(x))
  percentVar <- round(p$sdev^2/sum(p$sdev^2) * 100, 1)
  ggdat <- cbind(cdata, as.data.frame(p$x))
  p1 <- ggplot(ggdat, aes(PC1, PC2, color=Group, label=rownames(ggdat))) + 
    geom_point(size=3) + 
    ggrepel::geom_text_repel() +
    xlab(paste0("PC1: ", percentVar[1],"% variance")) + 
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
  ggsave(filename = outfile, plot = p1, width = 6, height = 6)
} 




# jaccard table
jaccardfile <- paste0(opt$input,'.jaccard.txt')
if(file.exists(jaccardfile)) {
    jdata <- read.table(jaccardfile,header=TRUE,sep="\t", row.names = 1, check.names = F)
    pheatmap(jdata, col=colorRampPalette(brewer.pal(9, 'Blues'))(255), 
             main = "Spatial Overlap Between Peaks", 
             filename = paste0(opt$input, '.jaccard.', opt$format), 
             width = 8, height = 8)
}



# number of intersections
intersectfile <- paste0(opt$input,'.n_intersection.txt')
if(file.exists(intersectfile)) {
    ndata <- read.table(intersectfile,header=TRUE,sep="\t", 
                        row.names = 1, check.names = F)
    pheatmap(ndata, col=colorRampPalette(brewer.pal(9, 'Greens'))(255), 
             main = "Number of Peak Intersections", 
             filename = paste0(opt$input, '.n_intersection.', opt$format), 
             width = 8, height = 8)
}



# qvalue table
qfile <- paste0(opt$input, "_qvalue.txt")
if(file.exists(qfile)) {
    qdata = read.table(qfile, header=TRUE, sep="\t", 
                       row.names = 1, check.names = F, na.strings = '.')
    if ( colnames(qdata)[1] == "Name" ) {
        # exclude this column name
        qdata <- qdata[,2:ncol(qdata)]
    }
    qdata[is.na(qdata)] <- 0
    makehm(qdata, 0, opt$qmax, colorRampPalette(brewer.pal(9, 'Reds'))(256), 
           "Mean Q-Value",
           paste0(opt$input,"_qvalue_hm"))
}



# log2 Fold Enrichment table
lfefile <- paste0(opt$input,"_log2FE.txt")
if(file.exists(lfefile)) {
    lfedata = read.table(lfefile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = F, na.strings = '.')
    if ( colnames(lfedata)[1] == "Name" ) {
        # exclude this column name
        lfedata <- lfedata[,2:ncol(lfedata)]
    }
    lfedata[is.na(lfedata)] <- 0
    
    clrs <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))(256)
    makehm(lfedata, opt$min, opt$max, clrs, 
           "Mean Log2 Fold Enrichment",
           paste0(opt$input,"_log2FE_hm"))
    makekmeans(lfedata, opt$min, opt$max, 6, clrs, color6, 
               "Mean Log2 Fold Enrichment, 6 Clusters",
               paste0(opt$input,"_log2FE_hm_K6"))
    makekmeans(lfedata, opt$min, opt$max, 8, clrs, color8, 
               "Mean Log2 Fold Enrichment, 8 Clusters",
               paste0(opt$input,"_log2FE_hm_K8"))
    makekmeans(lfedata, opt$min, opt$max, 10, clrs, color10, 
               "Mean Log2 Fold Enrichment, 10 Clusters",
               paste0(opt$input,"_log2FE_hm_K10"))
}



# count data
countfile <- paste0(opt$input,"_counts.txt")
samplefile <- paste0(opt$input,"_samples.txt")
if(file.exists(countfile)) {
    # read input data
    countdata <- read.table(countfile, header=TRUE, sep="\t", row.names = 1, 
                            check.names = FALSE, na.strings = ".")
    sampledata <- read.table(samplefile, header=TRUE, sep="\t", row.names = 1, 
                             check.names = FALSE)
    
    if ( colnames(countdata)[1] == "Name" ) {
        # exclude this column name
        countdata <- countdata[,2:ncol(countdata)]
    }
    countdata[is.na(countdata)] <- 0
    
    # plot
    makeDist(countdata,sampledata,paste0(opt$input,"_distance.", opt$format))
    makePearCorr(countdata,sampledata,paste0(opt$input,"_pearson.", opt$format))
    makeSpearCorr(countdata,sampledata,paste0(opt$input,"_spearman.", opt$format))
    makePCA(countdata,sampledata,paste0(opt$input,"_PCA.", opt$format))
}




