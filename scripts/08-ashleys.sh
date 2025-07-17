#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Run ASHLEYS automated QC on directory of bam files."
    echo "OUTPUT:  {bamdirname}.features.tsv {bamdirname}.libquality.txt {bamdirname}.features_window_distribution.tsv {bamdirname}.libquality.log"
    echo ""
    echo "USAGE: bash $0 -i <bam directory> -g <genome size>"
    echo -e "\t-i path       Input bam directory. *Required*"
    echo -e "\t-g path       Genome size."
    echo -e "\t-t int        Number of threads. [Default: 12]"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:g:t:h" opt
do
    case "$opt" in
        i ) bamdir="$OPTARG" ;;
        g ) genomesize="$OPTARG" ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $bam || -z $genomesize ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $threads ]]; then
    threads=12 
fi

## SCRIPT
sampleId=$(basename $bamdir)

ashleys.py -j $threads features -f $bamdir -w 5000000 2000000 1000000 800000 600000 400000 200000 -o $sampleId.features.tsv
ashleys.py predict -p $sampleId.features.tsv -o $sampleId.libquality.txt -m /projects/lansdorp/nextflow_pipelines/alignStrandSeq/scripts/files/svc_default.pkl