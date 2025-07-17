#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Calling SNPs from the bam directory."
    echo "OUTPUT:  ./{sample}.bcftools.dp{depth}.vcf.gz ./{sample}.bcftools.dp{depth}.vcf.gz.tbi"
    echo ""
    echo "USAGE: bash $0 -i <bam directory> -g <genome> -d <depth>"
    echo -e "\t-i path       Input bam directory. *Required*"
    echo -e "\t-g path       Reference genome. *Required*"
    echo -e "\t-d int        Minimum depth of reads for SNP calling. [Default: 3]"
    echo -e "\t-t int        Number of threads. [Default: 12]"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:g:d:t:h" opt
do
    case "$opt" in
        i ) bamdir="$OPTARG" ;;
        g ) genome="$OPTARG" ;;
        d ) depth="$OPTARG" ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $bamdir || -z $genome ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $depth ]]; then
    depth=3 
fi
if [[ -z $threads ]]; then
    threads=12 
fi

## SCRIPT
ref=$(ls /projects/lansdorp/references/$genome/*.fa)
sample=$(basename $bamdir)
vcffile=$sample.bcftools.dp$depth.vcf.gz

seq 1 22 | sed 's/^/chr/' | shuf | xargs -P 4 -I {} sh -c "bcftools mpileup -Q 40 -q 30 -f $ref -r {} $bamdir/*bam | bcftools call --skip-variants indels -mv | bcftools view -i 'DP >= $depth' -Oz -o tmp.{}.vcf.gz"
bcftools concat -Oz -o $vcffile tmp.chr*.vcf.gz 
bcftools index $vcffile

if [[ -s $vcffile ]]; then
    rm tmp.chr*.vcf.gz
fi