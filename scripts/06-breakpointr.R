#!/usr/bin/env Rscript

## PURPOSE: BreakpointR for a specific directory.
## USAGE:   Rscript 06-breakpointr.R <paired> <threads> <bamdir>
## OUTPUT:  ./{bamdirId}_BPR_output/

## VARIABLES
args <- commandArgs(T)
paired <- as.logical(toupper(args[1]))
threads <- as.numeric(args[2])
bamdir <- args[3]

outdir <- paste0(basename(bamdir), "_BPR_output")

chr <- paste0("chr", c(seq(1,22), "X", "Y"))
mask <- "/projects/lansdorp/nextflow_pipelines/alignStrandSeq/scripts/files/blacklist.highdepth.centromeres.bed"

## SCRIPT
suppressMessages(suppressWarnings(library(breakpointR)))

breakpointr(inputfolder=bamdir, outputfolder=outdir, pairedEndReads=paired, numCPU=threads, 
    windowsize=175, binMethod="reads", background=0.15, chromosomes=chr, maskRegions=mask
)