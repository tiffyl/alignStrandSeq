#! /bin/bash

## PURPOSE: Output alignment rate based on raw bam file.
## USAGE:   bash stats-alignrate.sh <bam>
## OUTPUT:  ${libId}.colalmet

## VARIABLES
bam=$1
# ref=$(ls /projects/lansdorp/references/$2/bowtie2/*.fa)

## SCRIPT
libId=$(basename $bam .trimmed.bam)

java -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics \
    -I $bam \
    -O $libId.colalmet.txt \
    --VALIDATION_STRINGENCY LENIENT
    # -R ${ref} \

