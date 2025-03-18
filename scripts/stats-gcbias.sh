#! /bin/bash

## PURPOSE: Output insert size based on processed bam file.
## USAGE:   bash stats-insertsize.sh <bam> <ref>
## OUTPUT:  ${libId}.gc_bias.txt ${libId}.gc_bias.pdf ${libId}.summary_metrics.txt

## VARIABLES
bam=$1
ref=$(ls /projects/lansdorp/references/$2/bowtie2/*.fa)

## SCRIPT
libId=$(basename $bam .processed.bam) 

java -jar /usr/picard/picard.jar CollectGcBiasMetrics \
    -I $bam \
    -O $libId.gc_bias.txt \
    -CHART $libId.gc_bias.pdf \
    -S $libId.gc_summary.txt \
    -R $ref \
    --VERBOSITY DEBUG
