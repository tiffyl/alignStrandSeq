#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Create chr1:169000000-248956422 background for a specific sample/directory."
    echo "OUTPUT:  ./{bamdirId}.background.bam/.bai ./{bamdirId}.directional.bam/.bai ./{bamdirId}.background.stats ./{bamdirId}.directional.stats"
    echo ""
    echo "USAGE: bash $0 -i <bam directory> -p <paired> -f <filechr1wwlist>"
    echo -e "\t-i path       Input bam directory. *Required*"
    echo -e "\t-p bool       Paired-end reads. *Required*"
    echo -e "\t-f path       File with a list of libraries where chr1 is WW. *Required*"
    echo -e "\t-t int        Number of threads. [Default: 12]"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:p:f:t:h" opt
do
    case "$opt" in
        i ) bam="$OPTARG" ;;
        p ) paired="$OPTARG" ;;
        f ) chr1wwlist="$OPTARG" ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $bamdir || -z $paired || -z $chr1wwlist ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $threads ]]; then
    threads=12 
fi

## SCRIPT
sampleId=$(basename $bamdir)

mkdir ./tmp/

# Extract reverse-oriented reads/pairs (background) and forward-oriented reads/pairs (directional reads---we randomly take 1% of them).
grep -f <(ls $bamdir/*.bam | sed 's@.*/@@g' | cut -f1 -d ".") $chr1wwlist | while read filename
do
    if $paired; then
        samtools view -bh -F4 -F1024 -f81 $bamdir/$filename.*.bam chr1:169000000-248956422 > ./tmp/$filename.directional.r1.tmp
		samtools view -bh -F4 -F1024 -f161 $bamdir/$filename.*.bam chr1:169000000-248956422 > ./tmp/$filename.directional.r2.tmp
        samtools view -bh -F4 -F1024 -f97 $bamdir/$filename.*.bam chr1:169000000-248956422 > ./tmp/$filename.background.r1.tmp
        samtools view -bh -F4 -F1024 -f145 $bamdir/$filename.*.bam chr1:169000000-248956422 > ./tmp/$filename.background.r2.tmp
    else
        samtools view -bh -F4 -F1024 -f16 -s 111.01 $bamdir/$filename.*.bam chr1:169000000-248956422 > ./tmp/$filename.directional.r1.tmp
        samtools view -bh -F4 -F1024 -F16 $bamdir/$filename.*.bam chr1:169000000-248956422 > ./tmp/$filename.background.r1.tmp
    fi
done

# Merging for overall background analysis
samtools merge -@ $threads $sampleId.background.bam ./tmp/*background.*.tmp
samtools merge -@ $threads $sampleId.directional.bam ./tmp/*directional.*.tmp
samtools index -@ $threads $sampleId.background.bam
samtools index -@ $threads $sampleId.directional.bam

echo -e "$(samtools view -c $sampleId.background.bam) $(($(samtools view -c $sampleId.directional.bam)*100))" > $sampleId.readcounts.txt