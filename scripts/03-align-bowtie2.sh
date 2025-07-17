#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Align bam files to given reference genome (converted & sort to bam) with bowtie2."
    echo "OUTPUT:  ./{libId}.trimmed.bam"
    echo ""
    echo "USAGE: bash $0 -1 <r1.fastq> -p <paired (true/false)>"
    echo -e "\t-1 path       R1 fastq. *Required*"
    echo -e "\t-2 path       R2 fastq."
    echo -e "\t-p bool       Paired-end reads. *Required*"
    echo -e "\t-g str        Reference genome. *Required*" 
    echo -e "\t-t int        Number of threads. [Default: 12]"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "1:2:p:g:t:h" opt
do
    case "$opt" in
        1 ) r1fq="$OPTARG" ;;
        2 ) r2fq="$OPTARG" ;;
        p ) paired="$OPTARG" ;;
        g ) genome="$OPTARG" ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $r1fq || -z $paired || -z $genome ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $threads ]]; then
    threads=12 
fi

## SCRIPT
libId=$(basename $r1fq | cut -f1 -d ".")
ref=$(ls /projects/lansdorp/references/$genome/bowtie2/*.fa)

if $paired;
then
    bowtie2 -x $ref -p $threads -1 $r1fq -2 $r2fq | samtools sort -@ $threads -o $libId.trimmed.bam - 
else
    bowtie2 -x $ref -p $threads -U $r1fq | samtools sort -@ $threads -o $libId.trimmed.bam - 
fi

samtools index $libId.trimmed.bam