#!/usr/bin/env Rscript
# script to plot bam2wig shift models

args <- commandArgs(TRUE)
out <- args[1]

library(ggplot2)
library(reshape2)


plotShiftModel <- function(gdata, outFile) {
  gdata_long <- melt(gdata, id = "Start", variable.name = "Strand")
  p <- ggplot(gdata_long, aes(Start, value, xlim(-450,450))) + 
    geom_line(aes(color = Strand)) + 
    scale_x_continuous(breaks = c(-400,-200,0,200,400), 
                       name = "Relative to peak") + 
    ylab("Mean Read Depth") +
    ggtitle("Shift model")
  ggsave(plot = p, filename = outFile, width = 6, height = 3)
}

plotCorrelations <- function(gdata, outFile) {
  o <- order(gdata[,2], decreasing = T)
  colnames(gdata)[2] <- "Pearson"
  p <- ggplot(gdata, aes(Shift, Pearson, xlim(0,500))) + 
    scale_x_continuous(breaks = c(0,100,200,300,400,500), 
                       name = "Shift value") + 
    geom_line() +
    ylab("Mean Pearson Correlation") +
    geom_vline(xintercept = gdata[o[1],1], color = c("red")) +
    ggtitle("Pearson Shift Correlations")
  ggsave(plot = p, filename = outFile, width = 4, height = 3)
}

# plot the model PDF plot
mdata <- read.table(paste0(out, "_model.txt"), header = T)
plotShiftModel(mdata, paste0(out, "_model.pdf"))

# plot the correlation PDF plot
cdata <- read.table(paste0(out, "_correlations.txt"), header = T)
plotCorrelations(cdata, paste0(out, "_correlations.pdf"))


