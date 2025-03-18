#! /bin/bash

## PURPOSE: Output insert size based on processed bam file.
## USAGE:   bash stats-insertsize.sh <bam>
## OUTPUT:  ${libId}.colinsert.txt ${libId}.insert_size.pdf

## VARIABLES
bam=$1

## SCRIPT
libId=$(basename $bam .processed.bam)

java -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
    -I $bam \
    -O $libId.colinsert.txt \
    -H $libId.insert_size.pdf \
    -M 0.5