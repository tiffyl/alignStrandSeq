#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Filter bam files with provided minimum insert size length."
    echo "OUTPUT:  {sample}.ins{size}.bam"
    echo ""
    echo "USAGE: bash $0 -i <bam> -s <insertsize>"
    echo -e "\t-i path       Input BAM file. *Required*"
    echo -e "\t-s int        Minimum insert size to filter. *Required*"
    echo -e "\t-t int        Number of threads. [Default: 12]"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:s:t:h" opt
do
    case "$opt" in
        i ) bam="$OPTARG" ;;
        s ) size="$OPTARG" ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $bam || -z $size ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $threads ]]; then
    threads=12 
fi

## SCRIPT
sample=$(basename $bam | cut -f1 -d ".")
filteredbam=$sample.ins$size.bam

sambamba view -h -F "template_length >= $size or template_length <= -$size" -t $threads -f bam $bam -o $filteredbam