#! /bin/bash

## PURPOSE: Create chr1:169000000-248956422 background for a specific sample/directory.
## USAGE:   bash 09-background.sh <paired> <threads> <bamdir> <library list> <reference>
## OUTPUT:  ./{bamdirId}.background.bam/.bai ./{bamdirId}.directional.bam/.bai ./{bamdirId}.background.stats ./{bamdirId}.directional.stats

## VARIABLES
paired=$1
threads=$2
bamdir=$3
chr1wwlist=$4

## SCRIPT
sampleId=$(basename $bamdir)

mkdir ./tmp/

# Extract reverse-oriented reads/pairs (background) and forward-oriented reads/pairs (directional reads---we randomly take 1% of them).
grep -f <(ls $bamdir/*.bam | sed 's@.*/@@g' | cut -f1 -d ".") $chr1wwlist | while read filename
do
    if $paired; then
        samtools view -bh -F4 -F1024 -f81 -s 111.01 $bamdir/$filename.*.bam chr1:169000000-248956422 > ./tmp/$filename.directional.r1.tmp
		samtools view -bh -F4 -F1024 -f161 -s 111.01 $bamdir/$filename.*.bam chr1:169000000-248956422 > ./tmp/$filename.directional.r2.tmp
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