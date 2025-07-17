#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Extract number of unique reads a from given directory."
    echo "OUTPUT:  {bamdirname}.uniqreads.tsv "
    echo ""
    echo "USAGE: bash $0 -i <bam directory> "
    echo -e "\t-i path       Input bam directory. *Required*"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:t:h" opt
do
    case "$opt" in
        i ) sampledir="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $sampledir ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

## SCRIPT
sampleId=$(basename $sampledir)

total=0

for bam in $sampledir/*".bam"
do
	## Counting reads excluding unmapped and duplicates
	count=$(samtools view -c -F4 -F1028 $bam)
	total=$(($total + $count))
done

printf "$sampleId\t$total\n" > "$sampleId.uniqreads.tsv"


