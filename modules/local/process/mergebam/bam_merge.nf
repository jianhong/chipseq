// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * Merge bam file for replicates
 */
process JO_MERGE_REP_BAM {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:name) }

    conda (params.conda ? "./environment.txt" : null)

    input:
    tuple val(name), val(input), val(single_end), path(bam), path(inputbam)
    val options

    output:
    tuple val(name), path("${name}.*.sorted.bam"), path("${name}.*.sorted.bam.bai"), emit: bam
    tuple val(name), path("*.bw"), emit: bw

    script:
    def ioptions         = initOptions(options)
    def singleExt        = (single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
    def extendReads      = single_end ? "${singleExt}" : '--extendReads'
    def pe               = single_end?"se":"pe"
    """
    samtools merge \\
        ${name}.${pe}.bam \\
        ${bam.join(' ')}
    samtools sort -o ${name}.${pe}.sorted.bam ${name}.${pe}.bam
    samtools index ${name}.${pe}.sorted.bam

    bamCoverage -b ${name}.${pe}.sorted.bam \\
       -o ${name}.${pe}.norm.CPM.bw \\
       --binSize 10  --normalizeUsing CPM ${extendReads}

    if [ "input" != "null"]; then
     samtools merge \\
        ${input}.bam \\
        ${inputbam.join(' ')}
     samtools sort -o ${input}.sorted.bam ${input}.bam
     samtools index ${input}.bam
     bamCompare -b1 ${name}.${pe}.sorted.bam \\
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
    bamCoverage -b ${name}.${pe}.sorted.bam \\
       -o ${name}.${pe}.norm.RPGC.bw \\
       --effectiveGenomeSize $params.deep_gsize \\
       --binSize 10  --normalizeUsing RPGC ${extendReads}
    fi 
    """
}