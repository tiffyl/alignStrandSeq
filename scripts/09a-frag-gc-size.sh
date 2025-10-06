#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Obtain GC content and size of reads in given bamfile."
    echo "OUTPUT: ./{bamId}.stats ./{bamId}.stats"
    echo ""
    echo "USAGE: bash $0 -i <bam> -p <paired> -g <genome>"
    echo -e "\t-i path       Input bam. *Required*"
    echo -e "\t-g str        Reference genome. *Required*"
    echo -e "\t-p bool       Paired-end reads. *Required*"
    echo -e "\t-t int        Number of threads. [Default: 12]"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:g:p:t:h" opt
do
    case "$opt" in
        i ) bam="$OPTARG" ;;
        g ) genome="$OPTARG" ;;
        p ) paired="$OPTARG" ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $bam || -z $genome || -z $paired ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $threads ]]; then
    threads=12 
fi

# Obtain GC content and size of reads
ref=$(ls /projects/lansdorp/references/$genome/bowtie2/*.fa)
stats=$(basename $bam .bam)".stats"

if $paired;
then
	#for overall background analysis
	samtools sort -@ $threads -n $bam | bedtools bamtobed -bedpe -i stdin | awk '$1==$4 && ($6-$2)*($6-$2)<1000*1000 {print $1,$2,$6}' OFS="\t" | sort -k1,1 -k2,2n -k3,3n | bedtools nuc -bed stdin -fi $ref | cut -f5,12 > $stats
else
    samtools sort -@ $threads -n $bam | bedtools bamtobed -i stdin | sort -k1,1 -k2,2n -k3,3n | bedtools nuc -bed stdin -fi $ref | cut -f5,12 > $stats
fi