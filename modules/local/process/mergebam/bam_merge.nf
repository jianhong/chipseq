// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * Merge bam file for replicates
 */
process JO_MERGE_REP_BAM {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:meta.id) }
    
    container "quay.io/biocontainers/deeptools:3.4.3--py_0"
    conda (params.conda ? "${params.modules_dir}/mergebam/environment.txt" : null)

    input:
    tuple val(meta), path(bam), path(inputbam)
    val options

    output:
    tuple val(meta), path("${meta.id}.*.sorted.bam"), path("${meta.id}.*.sorted.bam.bai"), emit: bam
    tuple val(meta), path("*.bw"), emit: bw

    script:
    def ioptions         = initOptions(options)
    def singleExt        = (meta.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
    def extendReads      = meta.single_end ? "${singleExt}" : '--extendReads'
    def pe               = meta.single_end?"se":"pe"
    """
    samtools merge \\
        ${meta.id}.${pe}.bam \\
        ${bam.join(' ')}
    samtools sort -o ${meta.id}.${pe}.sorted.bam ${meta.id}.${pe}.bam
    samtools index ${meta.id}.${pe}.sorted.bam

    bamCoverage -b ${meta.id}.${pe}.sorted.bam \\
       -o ${meta.id}.${pe}.norm.CPM.bw \\
       --binSize 10  --normalizeUsing CPM ${extendReads}

    if [ "${meta.control}" != "null"]; then
     samtools merge \\
        ${meta.control}.bam \\
        ${inputbam.join(' ')}
     samtools sort -o ${meta.control}.sorted.bam ${meta.control}.bam
     samtools index ${meta.control}.bam
     bamCompare -b1 ${meta.id}.${pe}.sorted.bam \\
                -b2 ${meta.control}.sorted.bam \\
                --scaleFactorsMethod readCount \\
                --operation log2 \\
                --normalizeUsing CPM \\
                --pseudocount 1 \\
                --skipZeroOverZero \\
                -o  ${meta.id}.normByInput.CPM.log2ratio.bw
    fi
        
    if [ "$params.deep_gsize" != "" ] && [ "$params.deep_gsize" != "false" ] && [ "$params.deep_gsize" != "null" ]
    then
    bamCoverage -b ${meta.id}.${pe}.sorted.bam \\
       -o ${meta.id}.${pe}.norm.RPGC.bw \\
       --effectiveGenomeSize $params.deep_gsize \\
       --binSize 10  --normalizeUsing RPGC ${extendReads}
    fi 
    """
}