#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Output insert size based on processed bam file."
    echo "OUTPUT:  ${libId}.colinsert.txt ${libId}.insert_size.pdf"
    echo ""
    echo "USAGE: bash $0 -i <bam>"
    echo -e "\t-i path       Input bam. *Required*"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:t:h" opt
do
    case "$opt" in
        i ) bam="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $bam ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

## SCRIPT
libId=$(basename $bam .processed.bam)

java -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
    -I $bam \
    -O $libId.colinsert.txt \
    -H $libId.insert_size.pdf \
    -M 0.5