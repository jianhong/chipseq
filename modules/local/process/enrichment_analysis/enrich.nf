// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

/*
 * enrichment analysis by clusterProfile and GSEA
 */
process JO_ENRICHMENTANALYSIS {
    tag "$meta.id"
    label 'process_long'
    label 'error_ignore'
    publishDir "${params.outdir}/${meta[0].peaktype}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.conda ? "${params.conda_softwares.rbase}" : null)

    input:
    path res
    val options

    output:
    path 'enrichment/*', emit: enrichment

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def ioptions  = initOptions(options)
    """
    install_packages.r ChIPpeakAnno clusterProfiler pathview biomaRt optparse
    ${workflow.projectDir}/modules/local/process/enrichment_analysis/enrich.r \\
        -s '${params.genome}' \\
        $ioptions.args
    """
}