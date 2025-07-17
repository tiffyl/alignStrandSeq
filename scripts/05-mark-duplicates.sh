#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Mark duplicates in bam file."
    echo "OUTPUT:  ./{libId}.processed.bam ./{libId}.processed.bam.bai ./{libId}.mdup.txt"
    echo ""
    echo "USAGE: bash $0 -i <bam>"
    echo -e "\t-i path       Input bam. *Required*"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:h" opt
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
libId=$(basename $bam | cut -f1 -d ".")

java -jar /usr/picard/picard.jar MarkDuplicates \
	-I ${bam} \
	-O ${libId}.processed.bam \
	-M ${libId}.mdup.txt \
	--VALIDATION_STRINGENCY LENIENT \
	--CREATE_INDEX true

# Change the naming of bam index 
mv ${libId}.processed.bai ${libId}.processed.bam.bai
