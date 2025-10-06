#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Output alignment rate based on raw bam file."
    echo "OUTPUT:  ${libId}.colalmet"
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
libId=$(basename $bam .trimmed.bam)

java -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics \
    -I $bam \
    -O $libId.colalmet.txt \
    --VALIDATION_STRINGENCY LENIENT

