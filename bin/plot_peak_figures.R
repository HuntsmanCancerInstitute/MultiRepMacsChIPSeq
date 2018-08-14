#!/usr/bin/env Rscript
# script for plotting collected peak information


# run with basename as first argument in target directory



args <- commandArgs(TRUE)
expname = args[1]


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



makekmeans <-function(rdata, hm_min, hm_max, k, hm_color, rowColor, outbase) {
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
  pheatmap(kdata[,c(1:n)], cluster_cols = F, color = hm_color, 
           border_color = NA, breaks = seq(hm_min,hm_max,length=256), 
           cluster_rows = F, show_rownames = F, 
           show_colnames = T, annotation_row = rowAnnot, 
           annotation_colors = list("Cluster" = rowColor),
           filename = paste0(outbase, ".png"), width = 8, height = 10)
  write.table(data.frame(kresult$cluster), file = paste0(outbase, "_klist.txt"), sep = "\t", quote = F)
}

makehm <-function(rdata, hm_min, hm_max, hm_color, outbase) {
  o <- apply(rdata, 1, mean)
  rdata <- rdata[order(o, decreasing = T),]
  rdata <- replace(rdata, rdata > hm_max, hm_max)
  rdata <- replace(rdata, rdata < hm_min, hm_min)
  pheatmap(rdata, cluster_cols = F, color = hm_color, 
           border_color = NA, breaks = seq(hm_min,hm_max,length=256), 
           cluster_rows = F, show_rownames = F, show_colnames = T, 
           annotation_colors = list("stage" = stage_color),
           filename = paste0(outbase, ".png"), width = 8, height = 10)
}


# jaccard table
jdata <- read.table(paste0(expname,'.jaccard.txt'),header=TRUE,sep="\t", row.names = 1, check.names = F)
pheatmap( jdata, col=colorRampPalette(brewer.pal(9, 'Blues'))(255), 
          main = "Jaccard metric of spatial overlap", 
          filename = paste0(expname, '.jaccard.pdf'), width = 8, height = 8)

# number of intersections
ndata <- read.table(paste0(expname,'.n_intersection.txt'),header=TRUE,sep="\t", row.names = 1, check.names = F)
pheatmap( ndata, col=colorRampPalette(brewer.pal(9, 'Greens'))(255), 
          main = "Number of intersections", 
          filename = paste0(expname, '.n_intersection.pdf'), width = 8, height = 8)

# qvalue table
qdata = read.table(paste0(expname,"_qvalue.txt"),header=TRUE,sep="\t", row.names = 1, check.names = F, na.strings = '.')
makehm(qdata, 0, 30, colorRampPalette(brewer.pal(9, 'Reds'))(256), paste0(expname,"_qvalue_hm"))

# log2FE table
lfedata = read.table(paste0(expname,"_log2FE.txt"),header=TRUE,sep="\t", row.names = 1, check.names = F, na.strings = '.')
makehm(lfedata, -4, 4, colorRampPalette(rev(brewer.pal(9, 'RdBu')))(256), paste0(expname,"_log2FE_hm"))
makekmeans(lfedata, -4, 4, 6, colorRampPalette(rev(brewer.pal(9, 'RdBu')))(256), color6, paste0(expname,"_log2FE_hm_K6"))
makekmeans(lfedata, -4, 4, 8, colorRampPalette(rev(brewer.pal(9, 'RdBu')))(256), color8, paste0(expname,"_log2FE_hm_K8"))
makekmeans(lfedata, -4, 4, 10, colorRampPalette(rev(brewer.pal(9, 'RdBu')))(256), color10, paste0(expname,"_log2FE_hm_K10"))




