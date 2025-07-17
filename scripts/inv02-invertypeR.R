#! /usr/bin/env Rscript

## PURPOSE: InvertypeR for a specific directory.
## USAGE:   Rscript inv02-invertyper.R -i <bamdir> -v <vcf> -p <paired>  
## OUTPUT:  ./{bamdirId}_invertypeR/

suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(invertyper)))

## VARIABLES
parser <- ArgumentParser()

parser$add_argument("-i", "--input_folder",
    type = "character", required = F, default=".",
    help = "Absolute path to the directory containing good-quality Strand-seq libraries for your sample. Default: '.'.",
    metavar = "/path/to/BAMs/"
)

parser$add_argument("-v", "--vcf",
    type = "character", 
    help = "Absolute path to a VCF file of SNVs to phase.",
    metavar = "/path/to/snvs.vcf"
)

parser$add_argument("-p", "--paired",
    type = "logical", required = F, default = TRUE,
    help = "Are the Strand-seq reads paired end? Default: TRUE.",
    metavar = "TRUE/FALSE"
)

parser$add_argument("-t", "--threads",
    type = "integer", required = F, default= 4,
    help = "The number of parallel threads to run on. Default: 4.",
    metavar = "INT"
)

parser$add_argument("--inversions",
    type = "character", required = F, 
    default=file.path("/projects/lansdorp/nextflow_pipelines/alignStrandSeq/scripts/files/hanlon_2021_inversions_BMCgenomics_augmented.bed"),
    help = "Absolute path to a BED file containing genomic intervals that might be inversions.",
    metavar = "/path/to/BED"
)

parser$add_argument("--hard_mask",
    type = "character", required = F, 
    default=file.path("/projects/lansdorp/nextflow_pipelines/alignStrandSeq/scripts/files/hard_mask_highdepth_centromeres.hg38.humans.bed"),
    help = "Absolute path to a BED file containing regions with unreliable Strand-seq data.",
    metavar = "/path/to/BED"
)

parser$add_argument("--soft_mask",
    type = "character", required = F, 
    default=file.path("/projects/lansdorp/nextflow_pipelines/alignStrandSeq/scripts/files/soft_mask.hg38.humans.bed"),
    help = "Absolute path to a BED file containing regions, like very large inversions, that
        occasionally interfere with composite file creation. Rarely really necessary (see InvertypeR documentation)",
    metavar = "/path/to/BED"
)

parser$add_argument("--prior",
    type = "character", required = F, 
    default = "0.9683,0.0158,0.0158", 
    help = "A comma-separated list of prior weights for inversion genotypes, without spaces. Only needs to be altered if --inversion_list is 
        not the default. E.g., 0.96,0.02,0.02. See InvertypeR for more details.",
    metavar = "comma-separated numbers"
)

args <- parser$parse_args()

chr <- paste0("chr", c(seq(1,22)))

samplename <- basename(args$input_folder)
outfoldername <- paste0(samplename, "_invertypeR")
args$prior <- as.numeric(unlist(strsplit(gsub(" ","",args$prior),",")))
args$paired <- as.logical(toupper(args$paired))

## SCRIPT
inversions <- invertyper_pipeline(
    regions_to_genotype = args$inversions,
    prior = args$prior,
    adjust_method = "all",
    input_folder = args$input_folder,
    output_folder = file.path("./", outfoldername),
    haploid_chromosomes = NULL,
    vcf = args$vcf,
    paired_reads = args$paired,
    confidence = 0.95,
    hard_mask = args$hard_mask,
    soft_mask = args$soft_mask,
    chromosomes = chr,
    numCPU = args$threads,
    save_composite_files = TRUE,
    write_browser_files = FALSE,
    discover_breakpointr_inversions = TRUE,
    breakpointr_prior = c(0.9, 0.05, 0.05),
    breakpointr_haploid_prior = c(0.9, 0.1),
    windowsize = c(120, 360),
    minReads = c(50, 50),
    background = 0.2,
    output_file = "inversions.txt"
)