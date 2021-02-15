// Import generic module functions
include { initOptions; saveFiles; getRealPath } from '../functions'

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
    
    val id
    val peaktype
    
    path ('macs2DiffBind/*')
    path ('homerDiffBind/*')
    
    path software_version
    val options

    output:
    path "index.html"
    
    script:
    def sampleLabel = id.join(' ')
    def pts         = peaktype.join(' ')
    def curr_path   = getRealPath()
    """
    n=($sampleLabel)
    s=($pts)
    paste <(printf '%s\n' "\${n[@]}") <(printf '%s\n' "\${s[@]}") > peaktype_files.txt
    cp ${index_docs} new.rmd
    ${curr_path}/utilities/install_packages.r rmarkdown DT ggplot2 reshape2
    Rscript -e "rmarkdown::render('new.rmd', output_file='index.html', params = list(design='${designtab}', genome='${params.genome}', species='${params.species}', summary='${workflow_summary}', exec_info='${params.tracedir}', launchdir='${workflow.launchDir}', conda='${params.conda}'))"
    """
}