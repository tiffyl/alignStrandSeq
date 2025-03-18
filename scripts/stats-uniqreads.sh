#! /bin/bash

## PURPOSE: Extract Number of good, bad libraries and unique reads a from given directory.
## USAGE:   bash stats-uniqreads.sh <bam directory> 
## OUTPUT:  {bamdirname}.uniqreads.tsv 

## VARIABLES
sampledir=$1

## SCRIPT
sampleId=$(basename $sampledir)

total=0

for bam in $sampledir/*".bam"
do
	## Counting reads excluding unmapped and duplicates
	count=$(samtools view -c -F4 -F1028 $bam)
	total=$(($total + $count))
done

printf "$sampleId\t$total\n" > "$sampleId.uniqreads.tsv"


