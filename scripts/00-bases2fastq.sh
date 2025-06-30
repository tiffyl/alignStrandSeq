#! /bin/bash

## PURPOSE: Run bases2fastq based on given run folder; move the fastqs into sample-specific directory in input folder of the run
## USAGE:   bash 00-bases2fastq.sh <threads> <rundir> (<mismatch threshold>)
## OUTPUT:  ./input ./b2fqc 

## VARIABLES
threads=$1
rundir=$(echo $2 | sed 's@/$@@g')
b2fdir=$rundir"_b2f"

## SCRIPT
if [[ -z $3 ]]; then
    bases2fastq -p $threads $rundir $b2fdir
else
    bases2fastq -p $threads --settings "I1MismatchThreshold,$3" --settings "I2MismatchThreshold,$3" $rundir $b2fdir 
fi

ls -d $b2fdir/Samples/* | grep -iv "Phix" | grep -iv "Unassigned" | while read sample; 
do
    mkdir -p ./input/$(basename $sample)
    mv $sample/*/*fastq.gz ./input/$(basename $sample)

    mkdir -p ./b2fqc
    mv $sample/*.html ./b2fqc 2>/dev/null
done

mkdir -p ./input/Unassigned
mv $b2fdir/Samples/Unassigned/*fastq.gz ./input/Unassigned

mv $b2fdir/info/Bases2Fastq.log ./b2fqc 
mv $b2fdir/*html ./b2fqc
