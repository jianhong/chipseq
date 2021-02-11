// Import generic module functions
include { saveFiles; getRealPath } from '../functions'

/*
 * Parse software version numbers
 */
process GET_SOFTWARE_VERSIONS {
    cache false
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:"pipeline_info", publish_id:'') }

    input:
    path versions
    val options

    output:
    path "software_versions.csv", emit: csv
    path 'software_versions_mqc.yaml', emit: yaml

    script:
    def curr_path = getRealPath()
    """
    echo $workflow.manifest.version > pipeline.version.txt
    echo $workflow.nextflow.version > nextflow.version.txt
    echo $workflow.manifest.name > pipelinename.version.txt
    echo $workflow.manifest.homePage > pipelienurl.version.txt
    
    ${curr_path}/get_software_versions/scrape_software_versions.py &> software_versions_mqc.yaml
    """
}
