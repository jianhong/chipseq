// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * Create index.html
 */
process JO_INDEX {
    tag "$name"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:name) }

    conda (params.conda ? "./environment.txt" : null)

    when:
    !params.skip_trackhub
    
    input:
    path index_docs
    path designtab
    path workflow_summary
    path images
    path doc_img
    path checksum
    path software_version
    val options

    output:
    path "trackhub/*"
    
    script:
    """
    cp ${index_docs} new.rmd
    Rscript -e "rmarkdown::render('new.rmd', output_file='index.html', params = list(peaktype='${PEAK_TYPE}', design='${designtab}', genome='${params.genome}', summary='${workflow_summary}'))"
    """
}