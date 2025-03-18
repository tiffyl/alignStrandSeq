#! /bin/bash

## PURPOSE: Concatenate fastqs from different lanes into a single fastq within provided directory
## USAGE:   bash 01-concat-lanes.sh <paired> <threads> <directory> 
## OUTPUT:  ./input/sample/*fastq.gz

## VARIABLES
paired=$1
threads=$2
dir=$3

## SCRIPT
for sample in $(ls $dir/*.fastq.gz | sed 's#.*/##g' | cut -f1 -d "_" | sed 's#-r..-c.*##g' | sort -u)
do
    mkdir -p ./input/$sample/

    ls $dir/$sample*.fastq.gz | sed 's#_L._R.*##g' | sort -u | parallel -j $threads cat {}*"L*_R1.fastq.gz" '>' ./input/$sample/{/}_R1.fastq.gz

    if $paired
    then
        ls $dir/$sample*.fastq.gz | sed 's#_L._R.*##g' | sort -u | parallel -j $threads cat {}*L*_R2.fastq.gz '>' ./input/$sample/{/}_R2.fastq.gz
    fi

done