#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Concatenate fastqs from different lanes into a single fastq within provided directory."
    echo "OUTPUT: ./input/sample/*fastq.gz"
    echo ""
    echo "USAGE: bash $0 -i <rundir> -p <paired (true/false)>"
    echo -e "\t-i path       Input Directory with fastqs. [Required]"
    echo -e "\t-p bool       Paired-end reads. [Required]"
    echo -e "\t-t int        Number of threads. [Default: 12]"
    exit 1 # Exit script after printing help
}

while getopts "i:p:t:" opt
do
    case "$opt" in
        i ) dir="$OPTARG" ;;
        p ) paired="$OPTARG" ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $dir ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $threads ]]; then
    threads=12 
fi

## SCRIPT
for sample in $(ls $dir/*.fastq.gz | sed 's#.*/##g' | cut -f1 -d "_" | sed 's#-r[[:digit:]]\+-c.*##g' | sort -u)
do
    mkdir -p ./input/$sample/

    ls $dir/$sample*.fastq.gz | sed 's#_L._R.*##g' | sort -u | parallel -j $threads cat {}*L*_R1.fastq.gz '>' ./input/$sample/{/}_R1.fastq.gz

    if $paired
    then
        ls $dir/$sample*.fastq.gz | sed 's#_L._R.*##g' | sort -u | parallel -j $threads cat {}*L*_R2.fastq.gz '>' ./input/$sample/{/}_R2.fastq.gz
    fi

done