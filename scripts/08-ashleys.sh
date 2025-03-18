#! /bin/bash

## PURPOSE: Run ASHLEYS automated QC on directory of bam files.
## USAGE:   bash 08-ashleys.sh <threads> <genomesize> <bam directory> 
## OUTPUT:  {bamdirname}.features.tsv {bamdirname}.libquality.txt {bamdirname}.features_window_distribution.tsv {bamdirname}.libquality.log 

## VARIABLES
threads=$1
genomesize=$2
bamdir=$3

## SCRIPT
sampleId=$(basename $bamdir)

ashleys.py -j $threads features -f $bamdir -w 5000000 2000000 1000000 800000 600000 400000 200000 -o $sampleId.features.tsv
ashleys.py predict -p $sampleId.features.tsv -o $sampleId.libquality.txt -m /projects/lansdorp/nextflow_pipelines/alignStrandSeq/scripts/files/svc_default.pkl