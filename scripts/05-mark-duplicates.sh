#! /bin/bash

## PURPOSE: Mark duplicates in bam file.
## USAGE:   bash 05-mark-duplicates.sh <bam>
## OUTPUT:  ./{libId}.processed.bam ./{libId}.processed.bam.bai ./{libId}.mdup.txt

## VARIABLES
bam=$1

## SCRIPT
libId=$(basename $bam .trimmed.filtered.bam)

java -jar /usr/picard/picard.jar MarkDuplicates \
	-I ${libId}.trimmed.filtered.bam \
	-O ${libId}.processed.bam \
	-M ${libId}.mdup.txt \
	--VALIDATION_STRINGENCY LENIENT \
	--CREATE_INDEX true

# Change the naming of bam index 
mv ${libId}.processed.bai ${libId}.processed.bam.bai
