#! /bin/bash

## PURPOSE: Filter bam files wiht provided minimum insert size length
## USAGE:   bash 04a-filterbam-insertsize.sh <threads> <bam> <insert size> 
##			Only for paired end reads as that's the only bam file that will provide insert size
## OUTPUT:  {sample}.ins{size}.bam

## VARIABLES
threads=$1
bam=$2
size=$3

## SCRIPT
sample=$(basename $bam | cut -f1 -d ".")
filteredbam=$sample.ins$size.bam

sambamba view -h -F "template_length >= $size or template_length <= -$size" -t $threads -f bam $bam -o $filteredbam