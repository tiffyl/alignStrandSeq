#! /bin/bash

## PURPOSE: Trim TruSeq adapters based on given R1 & R2 (optional) path
## USAGE:   bash 02-trim-adapters.sh <paired> <elevate> <r1.fastq> (<r2.fastq>) 
## OUTPUT:  ./{libId}_trimmed.1.fastq.gz (./{libId}_trimmed.2.fastq.gz) ./{libId}.trimming.log

## VARIABLES
paired=$1
elevate=$2
r1fq=$3
r2fq=$4

## SCRIPT
libId=$(basename $r1fq | cut -f1 -d "_")

if $elevate;
then
    r1_5adapt="CATGTAATGCACGTACTTTCAGGGTNNNNNNNNNCGTGCTGGATTGGCTCACCAGACACCTTCCGACAT"
    r1_3adapt="ATGTCGGAAGGTGTGCAGGCTACCGCTTGTCAACTNNNNNNNNNNNNAGTCGTCGCAGCCTCACCTGATC"

    r2_5adapt="GATCAGGTGAGGCTGCGACGACTNNNNNNNNNNNNAGTTGACAAGCGGTAGCCTGCACACCTTCCGACAT"
    r2_3adapt="ATGTCGGAAGGTGTCTGGTGAGCCAATCCAGCACGNNNNNNNNNACCCTGAAAGTACGTGCATTACATG"
else
    r1_5adapt="AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    r1_3adapt="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"

    r2_5adapt="CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
    r2_3adapt="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT"
fi

if $paired;
then
    cutadapt \
    -g $r1_5adapt \
    -a $r1_3adapt \
    -G $r2_5adapt \
    -A $r2_3adapt \
    -m 30 -q 15 \
    -o $libId".trimmed.1.fastq.gz" \
    -p $libId".trimmed.2.fastq.gz" \
    $r1fq $r2fq &> $libId".trimming.log"

else
    cutadapt \
    -g $r1_5adapt \
    -a $r1_3adapt \
    -m 30 -q 15 \
    -o $libId".trimmed.1.fastq.gz"  \
    $r1fq &> $libId".trimming.log"
fi
