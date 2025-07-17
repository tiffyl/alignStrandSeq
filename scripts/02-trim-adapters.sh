#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Trim TruSeq adapters based on given R1 & R2 (optional) path."
    echo "OUTPUT:  ./{libId}_trimmed.1.fastq.gz (./{libId}_trimmed.2.fastq.gz) ./{libId}.trimming.log"
    echo ""
    echo "USAGE: bash $0 -1 <r1.fastq> -p <paired (true/false)>"
    echo -e "\t-1 path       R1 fastq. *Required*"
    echo -e "\t-2 path       R2 fastq."
    echo -e "\t-p bool       Paired-end reads. *Required*"
    echo -e "\t-e            Elevate adapters trimming."
    echo -e "\t-n            Nextera adapters trimming."
    echo -e "\t-t int        Number of threads. [Default: 12]"
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "1:2:p:e:n:t:h" opt
do
    case "$opt" in
        1 ) r1fq="$OPTARG" ;;
        2 ) r2fq="$OPTARG" ;;
        p ) paired="$OPTARG" ;;
        e ) elevate="$OPTARG" ;;
        n ) nextera="$OPTARG" ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $r1fq || -z $paired ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $elevate ]]; then
    elevate=false  
fi
if [[ -z $nextera ]]; then
    nextera=false 
fi
if [[ -z $threads ]]; then
    threads=12 
fi

## SCRIPT
libId=$(basename $r1fq | cut -f1 -d "_")

if $elevate;
then
    r1_5adapt="CATGTAATGCACGTACTTTCAGGGTNNNNNNNNNCGTGCTGGATTGGCTCACCAGACACCTTCCGACAT"
    r1_3adapt="ATGTCGGAAGGTGTGCAGGCTACCGCTTGTCAACTNNNNNNNNNNNNAGTCGTCGCAGCCTCACCTGATC"

    r2_5adapt="GATCAGGTGAGGCTGCGACGACTNNNNNNNNNNNNAGTTGACAAGCGGTAGCCTGCACACCTTCCGACAT"
    r2_3adapt="ATGTCGGAAGGTGTCTGGTGAGCCAATCCAGCACGNNNNNNNNNACCCTGAAAGTACGTGCATTACATG"
elif $nextera;
then
    r1_5adapt="AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
    r1_3adapt="CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"

    r2_5adapt="CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    r2_3adapt="CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT"
else ## TruSeq
    r1_5adapt="AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    r1_3adapt="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"

    r2_5adapt="CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
    r2_3adapt="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT"
fi

if $paired;
then
    cutadapt \
    -g $r1_5adapt \
    -a $r1_3adapt \
    -G $r2_5adapt \
    -A $r2_3adapt \
    -m 30 -q 15 \
    -o $libId".trimmed.1.fastq.gz" \
    -p $libId".trimmed.2.fastq.gz" \
    $r1fq $r2fq &> $libId".trimming.log"

else
    cutadapt \
    -g $r1_5adapt \
    -a $r1_3adapt \
    -m 30 -q 15 \
    -o $libId".trimmed.1.fastq.gz"  \
    $r1fq &> $libId".trimming.log"
fi
