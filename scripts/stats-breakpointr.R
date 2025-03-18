#!/usr/bin/env Rscript

## PURPOSE: Obtain stats from BreakpointR RData.
## USAGE:   Rscript stats-breakpointr.R <BPR_directory>
## OUTPUT:  ./{sampleId}.bprstats.txt ./{sampleId}.wwchr1.list ./{sampleId}.chrstates.txt

## VARIABLES
args <- commandArgs(T)
bprdir <- args[1]

## SCRIPT
suppressMessages(suppressWarnings(library(breakpointR)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))

sampleId <- basename(strsplit(bprdir, "_BPR_output")[[1]][1])

filelist <- list.files(path=paste0(bprdir,"/data"), pattern="RData$", full.names=TRUE)
rdata <- loadFromFiles(filelist)

Library <- unname(unlist(lapply(rdata, function(x) strsplit(x$ID, "\\.")[[1]][1])))
Background <- unname(unlist(lapply(rdata, function(x) x$lib.metrics[1])))
Reads_per_Mb <- unname(unlist(lapply(rdata, function(x) x$lib.metrics[2])))
Percent_WC <- unname(unlist(lapply(rdata, function(x) length(grep("wc", x$counts$states)) / length(x$counts$states))))

stats <- data.frame(Library, Background, Reads_per_Mb, Percent_WC)
write.table(stats, paste0(sampleId, ".bprstats.txt"), quote=FALSE, row.names=FALSE, sep="\t")

# Output states of each chromosome
get_states <- function(bpr_object){
    df_obj <- as.data.frame(bpr_object$counts[width(bpr_object$counts) >= 10000000]) %>% 
        select(seqnames, states) %>% 
        group_by(seqnames) %>% 
        summarize(state = ifelse(n_distinct(states) == 1, first(states), NA)) %>% 
        pivot_wider(names_from = seqnames, values_from = state) %>%
        mutate(library = strsplit(bpr_object$ID, "\\.")[[1]][1]) %>%
        select(last_col(), sort(everything()))

    return(df_obj)
}

summary_states <- bind_rows(unname(lapply(rdata, get_states)))
write.table(summary_states, paste0(sampleId, ".chrstates.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

# Extract libraries with WW regions on chr1 (downstream background visualization)
bam_ww_chr1 <- summary_states %>% filter(chr1 == "ww") %>% select(library) %>% unlist() %>% unname() 
write.table(bam_ww_chr1, paste0(sampleId, ".wwchr1.list"), quote=FALSE, row.names=FALSE)

# get_ww_chr1 <- function(bpr_object){
#     test_region <- GRanges("chr1:169000000-248956422")
#     test_states <- subsetByOverlaps(bpr_object$counts, test_region)

#     if (length(test_states)==1 && test_states$states == "ww"){
#         return(strsplit(bpr_object$ID, "\\.")[[1]][1])
#     }
# }

# bam_ww_chr1 <- unname(unlist(lapply(rdata, get_ww_chr1)))
# write.table(bam_ww_chr1, paste0(sampleId, ".wwchr1.list"), quote=FALSE, row.names=FALSE)