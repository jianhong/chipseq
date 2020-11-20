// Import generic module functions
include { initOptions; saveFiles } from './functions'

/*
 * Merge bam file for replicates
 */
process MERGE_REP_BAM {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:name) }

    conda (params.conda ? "./environment.txt" : null)

    input:
    tuple val(name), val(input), path(bam), path(inputbam)
    val options

    output:
    tuple val(name), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    tuple val(name), path("*.bw"), emit: bw

    script:
    def ioptions         = initOptions(options)
    def singleExt        = (params.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
    def extendReads      = params.single_end ? "${singleExt}" : '--extendReads'
    """
    samtools merge \\
        ${name}.bam \\
        ${bam.join(' ')}
    samtools sort -o ${name}.sorted.bam ${name}.bam
    samtools index ${name}.sorted.bam

    bamCoverage -b ${name}.sorted.bam \\
       -o ${name}.norm.CPM.bw \\
       --binSize 10  --normalizeUsing CPM ${extendReads}

    if [ "input" != "null"]; then
     samtools merge \\
        ${input}.bam \\
        ${inputbam.join(' ')}
     samtools sort -o ${input}.sorted.bam ${input}.bam
     samtools index ${input}.bam
     bamCompare -b1 ${name}.sorted.bam \\
                -b2 ${input}.sorted.bam \\
                --scaleFactorsMethod readCount \\
                --operation log2 \\
                --normalizeUsing CPM \\
                --pseudocount 1 \\
                --skipZeroOverZero \\
                -o  ${name}.normByInput.CPM.log2ratio.bw
    fi
        
    if [ "$params.deep_gsize" != "" ] && [ "$params.deep_gsize" != "false" ] && [ "$params.deep_gsize" != "null" ]
    then
    bamCoverage -b ${name}.sorted.bam \\
       -o ${name}.norm.RPGC.bw \\
       --effectiveGenomeSize $params.deep_gsize \\
       --binSize 10  --normalizeUsing RPGC ${extendReads}
    fi 
    """
}