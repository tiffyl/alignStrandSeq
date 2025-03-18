# NEXTFLOW WORKFLOW FOR STRANDSEQ

Pipeline includes steps of bases2fastq, adapters trimming, alignment, stats collection.

<br>

## Running Workflow:
```
nextflow /projects/lansdorp/nextflow_pipelines/alignStrandSeq/wf_strandseq.nf --paired <true/false>
```

#### Additional Parameters:
```
--b2f <path>            Full path to run directory for bases2fastq with default parameters.
--b2fdir <path>         Full path to existing output directory of bases2fastq execution (Ran with default parameters).
--pipedir <path>        Full path to directory for pipeline (for input and output). [default: currect directory]
--ref <str>             Reference genome (Options in /projects/lansdorp/references/ with bowtie2 directory). [default: hg38]
--threads <int>         Number of threads. [default: 24]
--keeptrim              Keep trimmed fastqs. [default: false]
--keeprawbam            Keep raw bam files. [default: false]

-N <str>                Email to receive workflow notification.
-resume                 Resume nextflow execution from previous execution.
```
For more information: ```nextflow /projects/lansdorp/nextflow_pipelines/alignStrandSeq/wf_strandseq.nf --help```
<br><br>

### Recommendation for Running:
1. Run all nextflow workflows within ```/projects/lansdrop/lansdorp_scratch```.
2. Create a new directory following our typical naming convention. ```mkdir YYYY_MM_DD_SEQUENCER_HIGHorLOW_PE75or150```. <br>
    Can check ```/projects/lansdorp/bioinformatics_storage/``` for typical naming convention.
3. Move into the newly created directory. ```cd YYYY_MM_DD_SEQUENCER_HIGHorLOW_PE75or150```
4. Run nextflow workflow with parameters of choice. <br> 
    Would recommend running with -N <email> parameter as nextflow will receive emails for when workflow ends (completion or error).
5. If workflow stops at anytime for any reason. After identifying & fixing the error, rerun the workflow with ```-resume``` parameter to start from where the workflow stopped (and not from the beginning).
<br><br>

## Output:
- ```./input/``` Raw fastq files, separated into samples.
- ```./bam/``` Processed bam files, separated into samples.
- ```./analysis/``` All analysis reports, including heatmaps, breakplots, preseq and summary.
- ```./b2fqc/``` Log files from bases2fastq (if ran as part of workflow).
- Nextflow Report: execution report that provides information regarding resources and run time. 
- Nextflow Trace: information regarding individual each process/job that was ran within the workflow. 
- DAG: Workflow diagram
- Nextflow Timeline


