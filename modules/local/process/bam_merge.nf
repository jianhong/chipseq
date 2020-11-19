// Import generic module functions
include { initOptions; saveFiles } from './functions'

/*
 * Merge bam file for replicates
 */
process MERGE_REP_BAM {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:meta.id) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    tuple val(meta), path(bam), path(bai)
    val options

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.bw" emit: bw

    script:
    def ioptions         = initOptions(options)
    def singleExt        = (params.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
    def extendReads      = params.single_end ? "${singleExt}" : '--extendReads'
    """
    samtools merge \\
        ${meta}.bam \\
        ${bam.join(' ')}
        
    samtools sort -o ${meta}.sorted.bam ${meta}.bam 

    samtools index ${meta}.sorted.bam
    
    bamCoverage -b ${meta}.sorted.bam \\
       -o ${meta}.norm.CPM.bw \\
       --binSize 10  --normalizeUsing CPM ${extendReads}

    if [ "$params.deep_gsize" != "" ] && [ "$params.deep_gsize" != "false" ] && [ "$params.deep_gsize" != "null" ]
    then
    bamCoverage -b ${meta}.sorted.bam \\
       -o ${meta}.norm.RPGC.bw \\
       --effectiveGenomeSize $params.deep_gsize \\
       --binSize 10  --normalizeUsing RPGC ${extendReads}
    fi 
    """
}