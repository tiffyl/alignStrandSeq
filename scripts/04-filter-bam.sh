#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Filter bam files for unmapped reads & unwanted chromosomes."
    echo "OUTPUT:  ./{libId}.trimmed.filtered.bam"
    echo ""
    echo "USAGE: bash $0 -i <bam> -p <paired (true/false)> -m <mapq>"
    echo -e "\t-i path       Input BAM file. *Required*"
    echo -e "\t-p bool       Paired-end reads. *Required*"
    echo -e "\t-m int        Mapping Quality filter. *Required*"
    echo -e "\t-s bool       Input BAM needs to be sorted. [Default: true]"
    echo -e "\t-t int        Number of threads. [Default: 12]"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:p:m:s:t:h" opt
do
    case "$opt" in
        i ) bam="$OPTARG" ;;
        p ) paired="$OPTARG" ;;        
        m ) mapq="$OPTARG" ;;
        s ) sort="$OPTARG" ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $bam || -z $paired || -z $mapq ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $sort ]]; then
    sort=true 
fi
if [[ -z $threads ]]; then
    threads=12 
fi

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
    samtools view -F2052 -q $mapq -h $sorted chr{1..22} chrX chrY | grep -v -E '@SQ.*chrUn|@SQ.*random|@SQ.*chrEBV' | samtools sort -@ $threads -n - | 
        samtools fixmate -@ $threads -O bam - - | samtools sort -@ $threads - | samtools view -bh -f1 -o $libId.trimmed.filtered.bam
else
    samtools view -F2052 -q $mapq -h $sorted chr{1..22} chrX chrY | grep -v -E '@SQ.*chrUn|@SQ.*random|@SQ.*chrEBV' | 
        samtools view -bh -o $libId.trimmed.filtered.bam
fi
