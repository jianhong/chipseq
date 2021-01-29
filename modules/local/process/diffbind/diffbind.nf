// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

/*
 * Differential analysis by DiffBind and annotation by ChIPpeakAnno
 */
process JO_DIFFBIND {
    tag "$meta.id"
    label 'process_long'
    label 'error_ignore'
    publishDir "${params.outdir}/${meta[0].peaktype}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.conda ? "${params.conda_softwares.rbase} ${params.conda_softwares.libv8}" : null)

    input:
    tuple val(meta), path(bams), path(peaks)
    path gtf
    path blacklist
    val options

    output:
    path 'DiffBind/*', emit: res

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def ioptions  = initOptions(options)
    blacklist_params = params.blacklist ? "-b ${blacklist}" : '-b FALSE'
    def metadata = new groovy.json.JsonBuilder(meta).toPrettyString()
    """
    echo '${metadata}' > designtab.txt
    install_packages.r rjson DiffBind ChIPpeakAnno rtracklayer ggplot2 GenomicFeatures optparse
    diffbind.r -d 'designtab.txt' \\
        -f ${bams.join(',')} \\
        -p ${peaks.join(',')} \\
        -g ${gtf} \\
        ${blacklist_params} \\
        -s '${params.genome}' \\
        $ioptions.args
    """
}