// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * Create index.html
 */
process JO_INDEX {
    tag "$name"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:".", publish_id:name) }

    conda (params.conda ? "${params.conda_softwares.pandoc} ${params.conda_softwares.rbase}" : null)
    
    input:
    path index_docs
    path designtab
    path workflow_summary
    path images
    path doc_img
    path nsc
    path rsc
    path frip
    path preseq_log
    path flagstat
    path checksum
    path software_version
    val options

    output:
    path "index.html"
    
    script:
    """
    cp ${index_docs} new.rmd
    install_packages.r rmarkdown DT ggplot2 reshape2
    Rscript -e "rmarkdown::render('new.rmd', output_file='index.html', params = list(design='${designtab}', genome='${params.genome}', summary='${workflow_summary}', exec_info='${params.outdir}', launchdir='${workflow.launchDir}'))"
    """
}