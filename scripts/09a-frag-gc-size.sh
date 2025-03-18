#! /bin/bash

## PURPOSE: Obtain GC content and size of reads in given bamfile.
## USAGE:   bash 09a-frag-gc-size.sh <paired> <threads> <bam> <reference>
## OUTPUT:  ./{bamId}.stats ./{bamId}.stats

## VARIABLES
ppaired=$1
threads=$2
bam=$3
ref=$(ls /projects/lansdorp/references/$4/bowtie2/*.fa)

# Obtain GC content and size of reads
stats=$(basename $bam .bam)".stats"

if $paired;
then
	#for overall background analysis
	samtools sort -@ $threads -n $bam | bedtools bamtobed -bedpe -i stdin | awk '$1==$4 && ($6-$2)*($6-$2)<1000*1000 {print $1,$2,$6}' OFS="\t" | sort -k1,1 -k2,2n -k3,3n | bedtools nuc -bed stdin -fi $ref | cut -f5,12 > $stats
else
    samtools sort -@ $threads -n $bam | bedtools bamtobed -i stdin | sort -k1,1 -k2,2n -k3,3n | bedtools nuc -bed stdin -fi $ref | cut -f5,12 > $stats
fi