// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * Create trackhub
 */
process JO_TRACKHUB {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:'') }

    conda (params.conda ? "./environment.txt" : null)

    when:
    !params.skip_trackhub
    
    input:
    tuple val(name), path(bw)
    path designtab
    val options

    output:
    path "trackhub/*"
    
    script:
    def sampleLabel = name.join('___')
    def bws         = bw.join('___')
    """
    create_trackhub.py trackhub  $params.species $params.email $designtab '.mLb.clN__.mRp.clN' --path_prefix '../../../../'
    """
}