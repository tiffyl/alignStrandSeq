#! /bin/bash

## PURPOSE: Filter bam files for unmapped reads & unwanted chromosomes.
## USAGE:   bash 04-filter-bam.sh <paired> <threads> <bam> <sortbambycoordinates?(true/false)>
## OUTPUT:  ./{libId}.trimmed.filtered.bam

## VARIABLES
paired=$1
threads=$2
bam=$3
sort=$4

## SCRIPT
libId=$(basename $bam | cut -f1 -d ".")

if $sort;
then
    sorted=$libId.trimmed.sorted.bam

    samtools sort -@ $threads -o $sorted $bam
    samtools index $sorted
else
    sorted=$bam
fi

## TODO IF GENOME IS NOT HUMAN, THEN WE DO NOT FILTER FOR CHROMOSOMES
if $paired;
then
    #Remove unmapped, poor quality, unwanted chromosomes
    samtools view -F2052 -q10 -h $sorted chr{1..22} chrX chrY | grep -v -E '@SQ.*chrUn|@SQ.*random|@SQ.*chrEBV' | samtools sort -@ $threads -n - | 
        samtools fixmate -@ $threads -O bam - - | samtools sort -@ $threads - | samtools view -bh -f1 -o $libId.trimmed.filtered.bam
else
    samtools view -F2052 -q10 -h $sorted chr{1..22} chrX chrY | grep -v -E '@SQ.*chrUn|@SQ.*random|@SQ.*chrEBV' | 
        samtools view -bh -o $libId.trimmed.filtered.bam
fi
