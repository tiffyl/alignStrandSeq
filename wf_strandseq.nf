#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// HELP MESSAGE
if ( params.help ) {
    help = """wf_strandseq.nf:
             |
             |Optional arguments:
             |  --b2f <path>                Full path to run directory for bases2fastq with default parameters.
             |  --b2fdir <path>             Full path to existing output directory of bases2fastq execution (Ran with default parameters).
             |  --paired <true/false>       Paired end reads. [default: ${params.paired}]
             |  --elevate                   Elevate adapters are used. [default: ${params.elevate}]
             |  --nextera                   Elevate adapters are used. [default: ${params.elevate}]
             |  --ref <str>                 Reference genome (Options in /projects/lansdorp/references/ with bowtie2 directory). [default: ${params.ref}]
             |  --threads <int>             Number of threads. [default: ${params.threads}]
             |  --pipedir <path>            Full path to directory for pipeline (for input and output). [default: "./"]
             |  --invertyper                Run invertypeR on all files and create ideograms. [default: ${params.invertyper}]
             |  
             |  Filters:
             |  --ashleysthreshold <float>  Minimum threshold for Ashleys' QC filtering (0-1). [default: ${params.ashleysthreshold}]
             |  --wcthreshold <float>       Maximum threshold for WC fraction filtering (0-1). [default: ${params.wcthreshold}]
             |  --bgthreshold <float>       Maximum threshold for background filtering (0-1). [default: ${params.bgthreshold}]
             |  --mapq <int>                Mapping Quality threshold. [default: ${params.mapq}] 
             |  --insertsize <int>          Minimum length of insert size. [default: ${params.insertsize}]
             |  --depth <int>               Minimum depth of reads for SNP calling. [default: ${params.depth}]
             |  
             |  Outputs:
             |  --keeptrim                  Keep trimmed fastqs. [default: ${params.keeptrim}]
             |  --keeprawbam                Keep raw bam files. [default: ${params.keeprawbam}]
             |  --keepbpr                   Keep BPR_output directories. [default: ${params.keepbpr}]
             |
             |Optional Nextflow arugments:
             |  -N <str>                Email to receive workflow notification.
             |  -resume                 Resume nextflow execution from previous execution.
             | 
             |SetUp:
             |Naming convention: library identification is before '_' (_R1, _R2), and sample indentification is before '-r??-c??'.
             |By default, the pipeline will look for a 'input' directory in the current working directory for fastqs.
             |The fastqs should be separated into different sample directories (./input/{sample}/*fastq.gz)
             |If fastqs were separated by lanes, put them into a 'lanes' directory, and not separated into subdirectories. (./lanes/*fastq.gz).
             | |-input                                               |-lanes
             |   |- sample1                                             |- sample1-r1-c1_L1_R1.fastq.gz
             |      |- sample1-r1-c1_R1.fastq.gz                        |- sample1-r1-c1_L1_R2.fastq.gz
             |      |- sample1-r1-c1_R2.fastq.gz                        |- sample1-r1-c1_L2_R1.fastq.gz
             |      |- sample1-r1-c2_R1.fastq.gz                        |- sample1-r1-c1_L2_R2.fastq.gz
             |      |- sample1-r1-c2_R2.fastq.gz                        |- sample2-r2-c1_L1_R1.fastq.gz
             |   |- sample2                                             |- sample2-r2-c1_L1_R2.fastq.gz
             |      |- sample2-r2-c1_R1.fastq.gz                        |- sample2-r2-c1_L2_R1.fastq.gz
             |      |- sample2-r2-c1_R2.fastq.gz                        |- sample2-r2-c1_L2_R2.fastq.gz
             """.stripMargin()

    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

// PRCOESSES
include {
    bases2fastq
    concat_lanes
    trim_adapter
    bt2_align
    filter_insertsize
    remove_empty
    picard_mkdup
    breakpointr
    ashleysqc
    libraries_filter
    uniqreads
    preseq
    preseq_indiv
    preseq_group
    snpcalling
    invertyper
    ideogram
    picard_collectalign
    picard_collectinsertgc
    stats_adapter
    stats_empty
    stats_alignment
    stats_processed
    stats_insertsize
    stats_gcbias
    stats_breakpointr
    stats_ashleys
    stats_uniqreads
    stats_complexity
    stats_background
    plots_background
    metrics_summary
    output_trimmedfastq
    output_rawbams
    output_bpr
} from './strandseq_processes.nf'

// REFERENCES
if ( params.genomesize == null ){
    error("ERROR: Selected reference genome is not available for the alignStrandSeq pipeline.")
}

// SET UP
log.info """\
    STRANDSEQ PIPELINE
    ===================================
    started at          : ${workflow.start}
    ---
    output directory    : ${params.pipedir}
    reference genome    : ${params.ref}
    paired-end reads    : ${params.paired}
    threads             : ${params.threads}
    """
    .stripIndent(true)

workflow {
    //Bases2Fastq
    if ( params.b2fdir ) {
        //If bases2fastq previously ran, then can feed directory right into pipeline
        totrimfq_ch = Channel.fromPath("${params.b2fdir}/Samples/*/*/*fastq.gz").concat(Channel.fromPath("${params.b2fdir}/Samples/*/*fastq.gz"))
            .map{ file -> [ file.baseName.split("_")[0], file ] }
            .groupTuple()
            .filter{ id, files -> !(id =~ /(?i)^phix.*/) }
            .map{ id, files -> [ id, files.sort() ] }
    } 
    else if ( ! params.b2f ) {
        //If bases2fastq is skipped, then we must check if the fastqs are lane-split
        lanesdir = "${params.pipedir}/lanes/"
        if (file(lanesdir).isDirectory()) {
            println("Concatenating Fastqs from Lanes")
            concat_lanes(Channel.fromPath(lanesdir, checkIfExists: true))
            totrimfq_ch = concat_lanes.out.fastqs.flatMap{ it }
                .map{ file -> [ file.baseName.split("_")[0], file ] }
                .groupTuple()
                .map{ id, files -> [ id, files.sort() ] }
        }
        else {
            totrimfq_ch = Channel.fromPath("${params.pipedir}/input/*/*fastq.gz")
                .map{ file -> [ file.baseName.split("_")[0], file ] }
                .groupTuple()
                .map{ id, files -> [ id, files.sort() ] }
        }
    } 
    else {
        bases2fastq(Channel.fromPath(params.b2f))
        totrimfq_ch = bases2fastq.out.fastqs.flatMap{ it }
            .map{ file -> [ file.baseName.split("_")[0], file ] }
            .groupTuple()
            .filter{ id, files -> !(id =~ /(?i)^phix.*/) }
            .map{ id, files -> [ id, files.sort() ] }
    }

    //Trim
    trim_adapter(totrimfq_ch)
    
    //Align
    bt2_align(trim_adapter.out.trimmedFastqs)

    if ( params.paired && params.insertsize > 0 ) {
        processedBams_ch = filter_insertsize(bt2_align.out.bams)
    }
    else {
        processedBams_ch = bt2_align.out.bams
    }

    //BAM process
    remove_empty(processedBams_ch)
    picard_mkdup(remove_empty.out.toMkdupBams)

    sampleBams_ch = picard_mkdup.out.processedBams.groupTuple().map{ id, files -> [id, files.flatten()] }

    //BreakPointR
    breakpointr(sampleBams_ch)

    //Ashleys or QC
    ashleysqc(sampleBams_ch)
    tofilter_ch = ashleysqc.out.bamDirs.concat(breakpointr.out.bprstats)
        .concat(ashleysqc.out.ashleysQscores)
        .groupTuple(by:0)
    libraries_filter(tofilter_ch)
    uniqreads(libraries_filter.out.bamDirs)

    //PreSeq
    preseq(picard_mkdup.out.processedBams.filter{ id -> (id =~ /-r..-c../) })
    preseq_indiv(preseq.out.preseqFiles)
    preseq_group(preseq.out.preseqFiles.map{ id, file1, file2 -> [file1, file2]}.flatten().collect(), libraries_filter.out.poorLibs.collect())

    //InvertypeR
    if ( params.b2fdir ) {
        snpcalling(libraries_filter.out.bamDirs.filter{ id, dir -> !(id =~ /(?i)^un*/) })
        inv_ch = libraries_filter.out.bamDirs.concat(snpcalling.out.vcfs).groupTuple(by:0)
        invertyper(inv_ch)
        ideogram(invertyper.out.inversions)
    }
    
    //STATS
    stats_adapter(trim_adapter.out.trimLogs.collect())

    picard_collectalign(bt2_align.out.rawbams)
    stats_alignment(picard_collectalign.out.alignmets.collect())
    
    stats_empty(remove_empty.out.emptyBams.collect()) 

    stats_processed(picard_mkdup.out.mkdupLogs.collect())

    picard_collectinsertgc(picard_mkdup.out.processedBams)
    stats_insertsize(picard_collectinsertgc.out.insertmets.collect())
    stats_gcbias(picard_collectinsertgc.out.gcmets.collect())

    stats_breakpointr(breakpointr.out.bprstats.map{ id, file -> file }.flatten().collect(), breakpointr.out.chrstates.collect())
    stats_ashleys(ashleysqc.out.ashleysQscores.map{ id, file -> file }.flatten().collect())
    stats_uniqreads(uniqreads.out.uniqReads.collect())
    stats_complexity(preseq_indiv.out.complexity.collect())

    //BACKGROUND
    stats_background(libraries_filter.out.bamDirs.filter{ id, dir -> !(id =~ /(?i)^un*/) }.map{ id, dir -> dir }.collect(), 
                    breakpointr.out.wwchr1.filter{ file -> !(file.baseName =~ /(?i)^un*/) }.collect())
    plots_background(stats_background.out.bgFiles)

    //SUMMARY
    metrics_ch = stats_adapter.out.metrics.concat(
        stats_alignment.out.metrics,
        stats_processed.out.metrics,
        stats_insertsize.out.metrics,
        stats_gcbias.out.metrics,
        stats_breakpointr.out.metrics,
        stats_ashleys.out.metrics,
        stats_complexity.out.metrics,
        stats_uniqreads.out.metrics
    ).collect()
    metrics_summary(metrics_ch)

    //OUTPUT
    if ( params.keeptrim ) {
        output_trimmedfastq(trim_adapter.out.trimmedFastqs.map{ id, files -> files.flatten() }.collect())
    }
    if ( params.keeprawbam ) {
        output_rawbams(bt2_align.out.rawbams.collect())
    }
    if ( params.keepbpr ) {
        output_bpr(breakpointr.out.bprDirs.map{ id, file -> file }.collect())
    }

}

