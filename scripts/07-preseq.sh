#! /bin/bash

## PURPOSE: Generate complexity curves using preseq gc_extrap
## USAGE:   bash 07-preseq.sh <bam> <genomesize> (<failedlistfile>)
## OUTPUT:  ./{libId}.preseq.estimate ./{libId}.preseq.log (filefailedlist)

## VARIABLES
bam=$1
genomesize=$2
failedlistfile=$3

## SCRIPT
libId=$(basename $bam | cut -f1 -d ".")

# Covert BAM files to preseq's input format (Mapped Read)
bam2mr $bam > $libId.mr

# Obtain complexity estimate
if [ $genomesize -gt 300000000 ]; then
    preseq gc_extrap -v -e 3000000000 -o $libId.preseq.estimate $libId.mr &> $libId.preseq.log
elif [ $genomesize -gt 20000000 ]; then
    preseq gc_extrap -v -e 60000000 -s 1000000 -o $libId.preseq.estimate $libId.mr &> $libId.preseq.log
else
    preseq gc_extrap -v -e 6000000 -s 100000 -o $libId.preseq.estimate $libId.mr &> $libId.preseq.log
fi

rm $libId.mr

if [ ! -z $failedlistfile ] && [ ! -f "$libId.preseq.estimate" ];
then
    printf "$lidId\n" >> $failedlistfile
fi


## NOTES
###########################################################################################################################################
# Glossary
# Seqencing effort (measured in Mb or Gb) = total number of MAPPED bases after sequencing, including duplicates
# Breadth of coverage = proportion of the haploid genome that is covered at least one read, including sites with high depth 
#   (i.e. the algorithm is naive to alignment errors)

# Preseq
# The complexity curve is a function that predicts the number of unique reads that would be obtained if given number of reads were sequences
# i.e, if a library is sequenced with 100,000 reads, the complexity curve might predict 80,000 unique reads.

# One of the related publications makes the point that complexity curves for different libraries have different shapes and can cross
#   presumably because of different types and amount of amplification bias (e.g. by GC content) as compared with locus dropout 
#   (e.g.due to pre-amplification bead cleans).

# This means that the best way to compare two libraries is to plot their complexity curves side-by-side
#   but it also means that we can simplify by comparing their complexity curves at a fixed point (i.e. a fixed sequencing effort). 
#   This point of comparison should be chosen carefully.
#######################################################################################################################################
# April 30 2021
# A read clipping analysis shows that Preseq's complexity metric is independent of read length
# Tested David Porubsky's NA12878 libraries with 100bp PE reads: used GATK ClipReads to hard clip the reads to 75 bp. 
# The resulting Preseq plots were effectively identical.
# However, SE vs PE reads probably behave differently, because in PE reads the whole fragment counts towards predicted coverage
########################################################################################################################################
# Nov 21 2022
# A fragment clipping analysis shows that Preseq's complexity metric scales correctly with fragment size
# Comparison of a control set of 20 BAM files to the same files with 50 bp shorter fragments
# Removed 50 bp from the mapped coordinates, the read sequence, and the base qualities in the .mr file
# This reduced read length in single cells by 20%, and there was a roughly 20% reduction in Preseq complexity as well.