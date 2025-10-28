#!/usr/bin/env Rscript

## PURPOSE: BreakpointR for a specific directory.
## USAGE:   Rscript 06-breakpointr.R <paired> <threads> <bamdir> <ref>
## OUTPUT:  ./{bamdirId}_BPR_output/

## VARIABLES
args <- commandArgs(T)
paired <- as.logical(toupper(args[1]))
threads <- as.numeric(args[2])
bamdir <- args[3]
ref <- args[4]

outdir <- paste0(basename(bamdir), "_BPR_output")

if ( ref == "hg38" ) {
    chr <- paste0("chr", c(seq(1,22), "X", "Y"))
    mask <- "/projects/lansdorp/nextflow_pipelines/alignStrandSeq/scripts/files/hard_mask_highdepth_centromeres.hg38.humans.bed"
} else if ( ref == "chm13" ) {
    chr <- paste0("chr", c(seq(1,22), "X", "Y"))
    mask <- NULL
} else if ( ref == "mm39" ) {
    chr <- paste0("chr", c(seq(1,19), "X", "Y"))
    mask <- NULL
} else if ( ref == "canFam3" ) {
    chr <- paste0("chr", c(seq(1,38), "X"))
    mask <- NULL
} else if ( ref == "MesAur1" ) {
    chr <- paste0("KB", seq(708127,708150), ".1")
    mask <- NULL
} else if ( ref == "bGalGal1" ) {
    chr <- paste0("NC_0525", seq(32,72), ".1")
    mask <- NULL
}
    
## SCRIPT
suppressMessages(suppressWarnings(library(breakpointR)))

breakpointr(inputfolder=bamdir, outputfolder=outdir, pairedEndReads=paired, numCPU=threads, 
    windowsize=175, binMethod="reads", background=0.15, chromosomes=chr, maskRegions=mask
)