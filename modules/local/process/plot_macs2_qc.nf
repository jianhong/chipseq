// Import generic module functions
include { initOptions; saveFiles } from './functions'

/*
 * Aggregated QC plots for peaks
 */
process PLOT_MACS2_QC {
    label 'process_medium'
    publishDir "${params.outdir}/${peaktype}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:'') }

    conda (params.conda ? "${params.conda_softwares.rbase}" : null)

    input:
    tuple val(peaktype), path(peaks)
    val options

    output:
    path '*.txt', emit: txt
    path '*.pdf', emit: pdf

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def ioptions  = initOptions(options)
    """
    install_packages.r optparse ggplot2 reshape2 scales
    plot_macs2_qc.r \\
        -i ${peaks.join(',')} \\
        $ioptions.args
    """
}
