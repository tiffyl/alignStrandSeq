#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bases2fastq {
    container "${projectDir}/singularity/bases2fastq.sif"
    cpus "${params.totalCores / 2}"
    memory "60 GB"
    
    publishDir "${params.pipedir}/", mode: 'copy', pattern: '{input,b2fqc}'

    input:  
        path(b2fdir)

    output:
        path("input/*/*.fastq.gz"), emit: fastqs
        tuple path("input"), path("b2fqc")

    script:
    """
    bash ${projectDir}/scripts/00-bases2fastq.sh ${params.totalCores / 2} $b2fdir
    """
}

process concat_lanes {  
    publishDir "${params.pipedir}/", mode: 'copy', pattern: "input/*"

    input:
        path(lanesdir)

    output:
        path("input/*")
        path("input/*/*.fastq.gz"), emit: fastqs

    script:
    """
    bash ${projectDir}/scripts/01-concat-lanes.sh ${params.paired} ${params.subthreads} ${lanesdir}
    """
}

process trim_adapter {
    container "${projectDir}/singularity/cutadapt.sif"

    input:
        tuple val(libId), path(fastq)
    
    output:
        tuple val(libId), path("${libId}.trimmed.{1,2}.fastq.gz"), emit: trimmedFastqs
        path "${libId}.trimming.log", emit: trimLogs

    script:
    """
    if ${params.paired};
    then
        bash ${projectDir}/scripts/02-trim-adapters.sh ${params.paired} ${params.elevate} ${fastq[0]} ${fastq[1]} 
    else
        bash ${projectDir}/scripts/02-trim-adapters.sh ${params.paired} ${params.elevate} ${fastq}
    fi
    """
}

process bt2_align {
    container "${projectDir}/singularity/bowtie2_samtools_bedtools.sif"
    label "light_mem"
    label "high_cpus"

    input:
        tuple val(libId), path(trimmedfastq)
    
    output:
        tuple path("${libId}.trimmed.bam"), path("${libId}.trimmed.bam.bai"), emit: rawbams
        tuple val(libId), path("${libId}.trimmed.filtered.bam"), emit: bams
    
    script:
    """
    if ${params.paired};
    then
        bash ${projectDir}/scripts/03-align-bowtie2.sh ${params.paired} ${params.subthreads} ${params.ref} ${trimmedfastq[0]} ${trimmedfastq[1]} 
    else
        bash ${projectDir}/scripts/03-align-bowtie2.sh ${params.paired} ${params.subthreads} ${params.ref} ${trimmedfastq}
    fi

    bash ${projectDir}/scripts/04-filter-bam.sh ${params.paired} ${params.subthreads} ${libId}.trimmed.bam false

    """
}

process remove_empty {
    container "${projectDir}/singularity/bowtie2_samtools_bedtools.sif"

    input:
        tuple val(libId), path(bam)

    output:
        tuple val(libId), path(bam), optional:true, emit: toMkdupBams 
        path("empty/${bam}"), optional:true, emit: emptyBams

    script:
    """
    limit=5000

    reads=\$(samtools view -c ${bam})

    if [ \$reads -lt \$limit ]
    then
        mkdir -p ./empty
        mv ${bam} ./empty
    fi
    """
}

process picard_mkdup {
    container "${projectDir}/singularity/picard.sif"
    maxForks = 12
    label "heavy_mem"
    label "high_cpus"

    input:
        tuple val(libId), path(bam)
    
    output:
        tuple val(sampleId), path("${libId}.processed.*"), emit: processedBams
        path("${libId}.mdup.txt"), emit: mkdupLogs

    script:
    sampleId = libId.split(/-r\d{2}-c\d{2}/)[0]
    """
    bash ${projectDir}/scripts/05-mark-duplicates.sh ${bam}
    """
}

process breakpointr {
    container "${projectDir}/singularity/strandseq_Rtools.sif"
    publishDir "${params.pipedir}/analysis/breaksPlot/", mode: 'copy', pattern: "*breaksPlot.pdf"
    label "light_mem"

    input:
        tuple val(sampleId), path(bam)

    output:
        tuple val(sampleId), path("${sampleId}_BPR_output"), emit: bprDirs
        path("${sampleId}*.bprstats.txt"), emit: bprstats
        path("${sampleId}*.chrstates.txt"), emit: chrstates
        path("${sampleId}.wwchr1.list"), emit: wwchr1
        path("${sampleId}.breaksPlot.pdf")
    
    script:
    """
    mkdir -p ./${sampleId}
    mv *.bam* ./${sampleId}

    Rscript ${projectDir}/scripts/06-breakpointr.R ${params.paired} ${params.threads} ./${sampleId} 
    Rscript ${projectDir}/scripts/stats-breakpointr.R ${sampleId}_BPR_output

    mv ${sampleId}_BPR_output/plots/breaksPlot.pdf ${sampleId}.breaksPlot.pdf
    """
}

process preseq {
    container "${projectDir}/singularity/preseq.sif"

    input:
        tuple val(sampleId), path(bam)
    
    output:
        tuple val(sampleId), path("*.preseq.estimate"), path("*.preseq.log"), optional:true, emit: preseqFiles

    script:
    """
    bash ${projectDir}/scripts/07-preseq.sh ${bam[0]} ${params.genomesize}
    """
}

process preseq_indiv {
    container "${projectDir}/singularity/py310_viz.sif"
    publishDir "${params.pipedir}/analysis/preseq/sc/", mode: 'copy', pattern: "*.preseq.pdf"
    label "high_cpus"

    input:
        tuple val(sampleId), path(preseqEst), path(preseqLog)
    
    output:
        path("*.complexity.txt"), emit: complexity
        path("*.preseq.pdf")

    script:
    """
    python ${projectDir}/scripts/07a-preseq-indiv.py ${preseqEst} ${preseqLog} ${params.genomesize}
    """
}

process preseq_group {
    container "${projectDir}/singularity/py310_viz.sif"
    publishDir "${params.pipedir}/analysis/preseq/", mode: 'copy', pattern: "*.preseq.pdf"
    publishDir "${params.pipedir}/analysis/", mode: 'copy', pattern: "preseq_all.pdf"

    input:
        path(preseqEst)

    output:
        path("*.pdf")

    script:
    """
    python ${projectDir}/scripts/07b-preseq-group.py ${params.genomesize} "./" 
    """
}

process ashleysqc {
    container "${projectDir}/singularity/ashleysqc.sif"
    label "light_mem"
    label "high_cpus"

    publishDir "${params.pipedir}/bam", mode: 'copy', pattern: "${sampleId}"

    input:
        tuple val(sampleId), path(bams)

    output:
        path("${sampleId}/"), emit: bamDirs
        path("${sampleId}.libquality.txt"), emit: ashleysQscores

    script:
    """
    mkdir -p ./${sampleId}/poor_quality
    mv ${sampleId}*.processed.bam* ./${sampleId}

    bash ${projectDir}/scripts/08-ashleys.sh ${params.threads} ${params.genomesize} ./${sampleId}

    awk '(\$3 <= 0.5) {print \$1}' ${sampleId}.libquality.txt | while read bamfile; do mv ./${sampleId}/\$bamfile* ./${sampleId}/poor_quality; done
    """
}

process uniqreads {
    container "${projectDir}/singularity/bowtie2_samtools_bedtools.sif"

    input:
        path(bamdir)

    output:
        path("*.uniqreads.tsv"), emit: uniqReads

    script:
    """
    bash ${projectDir}/scripts/stats-uniqreads.sh ${bamdir}
    """
}

// OUTPUT FILES
process output_trimmedfastq {
    input:
        path(trimmedfastqs)

    script:
    """
    echo "Keeping Trimmed Fastqs" 
    
    mkdir -p ${params.pipedir}/trimmed_fastq
    cp *trimmed.*.fastq.gz ${params.pipedir}/trimmed_fastq
    """
}

process output_rawbams {
    input:
        path(rawbams)

    output:
        path("*trimmed.bam*")

    script:
    """
    echo "Keeping Raw Bams" 

    mkdir -p ${params.pipedir}/rawbams
    cp *trimmed.bam* ${params.pipedir}/rawbams
    """
}

process output_bpr {
    input:
        path(bprdir)

    output:
        path(bprdir)
    
    script:
    """
    echo "Keeping BPR Output Directories" 

    mkdir -p ${params.pipedir}/BPRs
    cp -Lr *BPR_output ${params.pipedir}/BPRs/
    """
}

// STATS PROCESSES
process picard_collectalign {
    container "${projectDir}/singularity/picard.sif"
    label "light_mem"
    label "high_cpus"

    input:
        tuple path(bam), path(bamidx)
    
    output:
        path("*colalmet.txt"), emit: alignmets

    script:
    """
    bash ${projectDir}/scripts/stats-alignment.sh ${bam}
    """    
}

process picard_collectinsertgc {
    container "${projectDir}/singularity/picard.sif"
    label "light_mem"
    label "high_cpus"

    publishDir "${params.pipedir}/analysis/insert_size/", mode: 'copy', pattern: "*.insert_size.pdf"
    publishDir "${params.pipedir}/analysis/gc_content/", mode: 'copy', pattern: "*.gc_bias.pdf"

    input:
        tuple val(sampleId), path(bam)
    
    output:
        path('*.colinsert.txt'), emit: insertmets
        path('*.insert_size.pdf')
        path('*.gc_bias.txt'), emit: gcmets
        path("*.gc_bias.pdf")
    
    """
    bash ${projectDir}/scripts/stats-insertsize.sh ${bam[0]} 
    bash ${projectDir}/scripts/stats-gcbias.sh ${bam[0]} ${params.ref}
    """
}

process stats_adapter {
    input:
        val(trimlogList)

    output:
        path("metrics_dimers.txt"), emit: metrics

    script:
    """
    printf "Library\tAdapters_at_least_50bp\n" > metrics_dimers.txt

    echo "${trimlogList.join("\n")}" | while read file
    do
        lib=\$(basename \$file .trimming.log)
        reads=\$(grep 'Total read' \$file | rev | cut -f1 -d" " | rev | sed 's/,//g')

        if [[ \$reads -gt 0 ]]
        then
            #Reporting the number of dimers found (>60bp)
            grep -E '^[0-9]{2}' \$file | awk -v reads="\$reads" -v lib=\$lib 'BEGIN {OFS="\t"} (\$1>=60){sum+=\$2}END{print lib,sum/(2*reads)}' >> metrics_dimers.txt
        fi
    done
    """
}

process stats_empty {
    publishDir "${params.pipedir}/analysis/", mode: 'copy', pattern: "empty.txt"

    input:
        val(emptyBamList)
    
    output:
        path('empty.txt') 
    
    script:
    """
        echo "${emptyBamList.join("\n")}" | sed 's@.*/@@g' | sed 's/\\..*//g' > empty.txt
    """
}

process stats_alignment {
    input:
        val(alignmetsList) 
    
    output:
        path("metrics_alignment.txt"), emit: metrics
    
    script:
    """
    printf "Library\tInitial_reads_aligned\tAlignment_rate\tRead_length\n" > metrics_alignment.txt
    
    if ${params.paired};
    then
        echo "${alignmetsList.join("\n")}" | while read file
        do
            lidId=\$(basename \$file .colalmet.txt)
            printf "\$lidId\t\$(grep CATEGORY -A3 \$file | tail -1 | cut -f6,7,16)\n" >> metrics_alignment.txt
        done
    else
        echo "${alignmetsList.join("\n")}" | while read file
        do
            lidId=\$(basename \$file .colalmet.txt)
            printf "\$lidId\t\$(grep CATEGORY -A1 \$file | tail -1 | cut -f6,7,16)\n" >> metrics_alignment.txt
        done
    fi
    """
}

process stats_processed {
    input:
        val(mkdupList) 
    
    output:
        path("metrics_processedstats.txt"), emit: metrics
    
    script:
    """
    printf "Library\tProcessed_reads_aligned\tDuplication_rate\n" > metrics_processedstats.txt

    if ${params.paired};
    then
        echo "${mkdupList.join("\n")}" | while read file
        do
            lidId=\$(basename \$file .mdup.txt)
            printf "\$lidId\t\$(grep LIBRARY -A1 \$file | tail -1 | cut -f3,9 | awk '{OFS="\t"}; {print \$1*2,\$2}')\n" >> metrics_processedstats.txt
        done
    else
        echo "${mkdupList.join("\n")}" | while read file 
        do
            lidId=\$(basename \$file .mdup.txt)
            printf "\$lidId\t\$(grep LIBRARY -A1 \$file | tail -1 | cut -f2,9)\n" >> metrics_processedstats.txt
        done
    fi
    """
}

process stats_insertsize {
    input:
        val(insertmetsList) 
    
    output:
        path("metrics_insertsize.txt"), emit: metrics
    
    script:
    """
    printf "Library\tMedian_insert_size\n" > metrics_insertsize.txt
    
    echo "${insertmetsList.join("\n")}" | while read file
    do
        libId=\$(basename \$file .colinsert.txt)
        printf "\$libId\t\$(grep -A1 "MEDIAN_INSERT_SIZE" \$file | tail -1 | cut -f1)\n" >> metrics_insertsize.txt 
    done
    """
}

process stats_gcbias {
    input:
        val(gcmetsList) 
    
    output:
        path("metrics_gc.txt"), emit: metrics
    
    script:
    """
    printf "Library\tMean_GC\n" > metrics_gc.txt
    
    echo "${gcmetsList.join("\n")}" | while read file
    do
        libId=\$(basename \$file .gc_bias.txt)
        printf "\$libId\t\$(grep "^All" \$file | cut -f3,5 | awk '{sum+=\$1*\$2; total+=\$2}; END {print sum/total/100}')\n" >> metrics_gc.txt 
    done
    """
}

process stats_breakpointr {
    publishDir "${params.pipedir}/analysis/", mode: 'copy', pattern: "libchromstates.txt"
    label "light_mem"
    
    input:
        val(bprstatsList)
        val(chromstatesList)
    
    output:
        path("metrics_bprstats.txt"), emit: metrics
        path("libchromstates.txt")
    
    script:
    """
    printf "Library\tBackground\tReads_per_Mb\tPercent_WC\n" > metrics_bprstats.txt
    cat ${bprstatsList.join(" ")} | grep -v "^Library" >> metrics_bprstats.txt

    printf "Library\t" > libchromstates.txt
    printf "chr%d\t" {1..22} >> libchromstates.txt
    printf "chrX\tchrY\n" >> libchromstates.txt
    cat ${chromstatesList.join(" ")} | sort >> libchromstates.txt
    """
}

process stats_ashleys {
    publishDir "${params.pipedir}/analysis/", mode: 'copy', pattern: "libquality.txt"

    input:
        val(ashleysList)
    
    output:
        path("metrics_libquality.txt"), emit: metrics
        path("libquality.txt")

    script:
    """
    printf "Library\tQuality\n" > metrics_libquality.txt
    cat ${ashleysList.join(" ")} | cut -f1,3 | grep -v "^cell" | sed 's/.processed.bam//g' >> metrics_libquality.txt

    cat ${ashleysList.join(" ")} | grep -v "^cell" | sort > libquality.txt
    """
}

process stats_uniqreads {
    publishDir "${params.pipedir}/analysis/", mode: 'copy'
    
    input:
        val(uniqReadsList)
    
    output:
        path("sample_uniqreads.tsv"), emit: metrics
    
    script:
    """
    printf "Sample\tGood_UniqReads\n" > sample_uniqreads.tsv
    cat ${uniqReadsList.join(" ")} | sort >> sample_uniqreads.tsv
    """
}

process stats_complexity {
    input:
        val(complexityList)
    
    output:
        path("metrics_complexity1Gb.txt"), emit: metrics
    
    script:
    """
    printf "Library\tComplexity_at_1Gb\n" > metrics_complexity1Gb.txt
    cat ${complexityList.join(" ")} >> metrics_complexity1Gb.txt
    """
}

process stats_background {
    container "${projectDir}/singularity/bowtie2_samtools_bedtools.sif"

    input:
        path(bamDirsList)
        path(wwchr1List)
    
    output:
        tuple path("*.stats"), path("*readcounts.txt"), emit: bgFiles

    script:
    """
    ls -d */ | while read bamdir
    do
        if ls \$bamdir/*bam* 1> /dev/null 2>&1; then
            bash ${projectDir}/scripts/09-background.sh ${params.paired} ${params.threads} \$bamdir \$(basename \$bamdir).wwchr1.list 
        fi
    done

    ## FOR ALL BACKGROUND
    samtools merge -@ ${params.threads} background.bam *background.bam 
    samtools merge -@ ${params.threads} directional.bam *directional.bam 
    samtools index -@ ${params.threads} background.bam
    samtools index -@ ${params.threads} directional.bam

    echo -e "\$(samtools view -c background.bam) \$((\$(samtools view -c directional.bam)*100))" > readcounts.txt

    
    ## Stats
    ls *.bam | while read bamfile;
    do
        bash ${projectDir}/scripts/09a-frag-gc-size.sh ${params.paired} ${params.threads} \$bamfile ${params.ref}
    done
    """
}

process plots_background {
    container "${projectDir}/singularity/py310_viz.sif"

    publishDir "${params.pipedir}/analysis/", mode: 'copy'

    input:
        tuple path(stats), path(readscounts)
    
    output:
        path("background.pdf")

    script:
    """
    python ${projectDir}/scripts/09b-background-plot.py ${params.paired}
    """    
}

process metrics_summary {
    container "${projectDir}/singularity/py310_viz.sif"
    publishDir "${params.pipedir}/analysis/", mode: 'copy'

    input:
        path(metricsFiles)

    output:
        path("metrics_details.tsv")
        path("heatmaps.pdf")
        path("metrics.xlsx")

    script:
    """
    python ${projectDir}/scripts/10-metrics_summary.py ${params.genomesize}
    python ${projectDir}/scripts/10a-metrics_plots.py
    """
}