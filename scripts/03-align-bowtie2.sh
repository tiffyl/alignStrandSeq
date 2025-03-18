#! /bin/bash

## PURPOSE: Align bam files to given reference genome (converted & sort to bam) with bowtie2
## USAGE:   bash 03-align-bowtie2.sh <paired> <threads> <reference> <r1.fastq> (<r2.fastq>) 
## OUTPUT:  ./{libId}.trimmed.bam

## VARIABLES
paired=$1
threads=$2
ref=$(ls /projects/lansdorp/references/$3/bowtie2/*.fa)
r1fq=$4
r2fq=$5

## SCRIPT
libId=$(basename $r1fq | cut -f1 -d ".")

if $paired;
then
    bowtie2 -x $ref -p $threads -1 $r1fq -2 $r2fq | samtools sort -@ $threads -o $libId.trimmed.bam - 
else
    bowtie2 -x $ref -p $threads -U $r1fq | samtools sort -@ $threads -o $libId.trimmed.bam - 
fi

samtools index $libId.trimmed.bam