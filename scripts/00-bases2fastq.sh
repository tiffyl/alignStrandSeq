#! /bin/bash

## PARAMETERS
helpFunction()
{
    echo "PURPOSE: Run bases2fastq based on given run folder; move the fastqs into sample-specific directory in input folder of the run."
    echo "OUTPUT:  ./input ./b2fqc"
    echo ""
    echo "USAGE: bash $0 -i <rundir>"
    echo -e "\t-i path       Input directory to run bases2fastq. *Required*"
    echo -e "\t-t int        Number of threads. [Default: 24]"
    echo -e "\t-m int        Number of mismatch for index demultiplexing. [Default: 1]"
    echo -e "\t-l            Split lanes for demutliplexing."
    echo -e "\t-h            Help message."

    exit 1 # Exit script after printing help
}

while getopts "i:mlt:h" opt
do
    case "$opt" in
        i ) rundir=$(echo $OPTARG | sed 's@/$@@g') ;;
        m ) mismatch="$OPTARG" ;;
        l ) splitlanes=1 ;;
        t ) threads="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Required
if [[ -z $rundir ]]; then
    echo "ERROR: Missing required parameters."
    helpFunction
fi

# Default
if [[ -z $threads ]]; then
    threads=24 
fi

## SCRIPT
b2fdir=$rundir"_b2f"

if [[ ! -z $mismatch && $mismatch != 1 && ! -z $splitlanes ]]; then
    bases2fastq -p $threads --settings "I1MismatchThreshold,$mismatch" --settings "I2MismatchThreshold,$mismatch" --split-lanes $rundir $b2fdir 
elif [[ ! -z $mismatch && $mismatch != 1 ]]; then
    bases2fastq -p $threads --settings "I1MismatchThreshold,$mismatch" --settings "I2MismatchThreshold,$mismatch" $rundir $b2fdir 
elif [[ ! -z $splitlanes ]]; then
    bases2fastq -p $threads --split-lanes $rundir $b2fdir 
else
    bases2fastq -p $threads $rundir $b2fdir
fi

ls -d $b2fdir/Samples/* | grep -iv "Phix" | grep -iv "Unassigned" | while read sample; 
do
    if [[ ! -z $splitlanes ]]; then
        mkdir -p ./L1/input/$(basename $sample) ./L2/input/$(basename $sample) 
        mv $sample/*/*_L1*fastq.gz ./L1/input/$(basename $sample)
        mv $sample/*/*_L2*fastq.gz ./L2/input/$(basename $sample)
    else
        mkdir -p ./input/$(basename $sample)
        mv $sample/*/*fastq.gz ./input/$(basename $sample)
    fi

    mkdir -p ./b2fqc
    mv $sample/*.html ./b2fqc 2>/dev/null
done

if [[ ! -z $splitlanes ]]; then
    mkdir -p ./L1/input/Unassigned ./L2/input/Unassigned
    mv $b2fdir/Samples/Unassigned/*_L1*fastq.gz ./L1/input/Unassigned
    mv $b2fdir/Samples/Unassigned/*_L2*fastq.gz ./L2/input/Unassigned
else
    mkdir -p ./input/Unassigned
    mv $b2fdir/Samples/Unassigned/*fastq.gz ./input/Unassigned
fi

mv $b2fdir/info/Bases2Fastq.log ./b2fqc 
mv $b2fdir/*html ./b2fqc
