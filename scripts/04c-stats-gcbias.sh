#! /bin/bash

## PURPOSE: Output insert size based on processed bam file.
## USAGE:   bash stats-insertsize.sh <bam> <ref>
## OUTPUT:  ${libId}.gc_bias.txt ${libId}.gc_bias.pdf ${libId}.summary_metrics.txt

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Output insert size based on processed bam file."
    echo "OUTPUT:  {libId}.gc_bias.txt {libId}.gc_bias.pdf {libId}.summary_metrics.txt"
    echo ""
    echo "USAGE: bash $0 -i <bam> -g <genome>"
    echo -e "\t-i path       Input bam. *Required*"
    echo -e "\t-g path       Reference genome. *Required*"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:g:t:h" opt
do
    case "$opt" in
        i ) bam="$OPTARG" ;;
        g ) genome="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $bam || -z $genome ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

## SCRIPT
ref=$(ls /projects/lansdorp/references/$genome/bowtie2/*.fa)
libId=$(basename $bam .processed.bam) 

java -jar /usr/picard/picard.jar CollectGcBiasMetrics \
    -I $bam \
    -O $libId.gc_bias.txt \
    -CHART $libId.gc_bias.pdf \
    -S $libId.gc_summary.txt \
    -R $ref \
    --VERBOSITY DEBUG
