#!/usr/bin/env Rscript

## PURPOSE: InvertypeR for a specific directory.
## USAGE:   Rscript inv03-ideogram.R <inversions> 
## OUTPUT:  ./inversions.png/

## VARIABLES
args <- commandArgs(T)
invfile <- args[1]

sample <- paste0(strsplit(invfile, '\\.')[[1]][1], ".inversions")
outfile <- paste0(sample, ".svg")

## SCRIPT
suppressMessages(suppressWarnings(library(RIdeogram)))

data(human_karyotype, package="RIdeogram")
human_karyotype <- human_karyotype[1:22,]
inversions <- read.table(invfile, sep = "\t", header=FALSE, col.names=c('Chr', 'Start', 'End', 'Value', 'Type', 'Shape'), stringsAsFactors = FALSE)

types <- inversions[,c('Type', 'Shape', 'Chr', 'Start', 'End')]
types$color <- 252525

ideogram(karyotype = human_karyotype, overlaid = inversions, label = types, label_type = "marker", output = outfile)
convertSVG(outfile, file = sample, device = "png")