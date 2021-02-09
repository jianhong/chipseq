// Import generic module functions
include { initOptions; saveFiles } from '../../functions'
params.options = [:]
/*
 * Create index.html
 */
process JO_SRAINFO {
    tag "$id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:".", publish_id:id) }

    conda (params.conda ? "${params.conda_softwares.entrez_direct}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "jianhong/chipseq:latest"
    }
    input:
    val id

    output:
    path "sampleInfo.csv", emit: info
    path "designTabFull.csv", emit: design
    path "*.version.txt", emit: version
    
    script:
    //def api_key = ${params.api_key} ? "-api_key ${params.api_key}": ""
    def api_key = ""
    """
    if [[ "${id}" == GSE* ]];then
        esearch -db gds -q "${id}" $api_key \\
            | elink -target sra $api_key \\
            | efetch -format runinfo $api_key > SraRunTable.${id}.txt
        esearch -query "${id}" -db gds $api_key \\
            | elink -target sra $api_key \\
            | efetch --json $api_key > SraRunInfo.${id}.txt
    else
        esearch -db sra -q "${id}" $api_key \\
            | efetch -format runinfo $api_key > SraRunTable.${id}.txt
        esearch -query "${id}" -db sra $api_key \\
            | efetch --json $api_key > SraRunInfo.${id}.txt
    fi
    
    ${workflow.projectDir}/modules/local/process/sra_downloader/sra_info/sra_samplesheet.py SraRunTable.${id}.txt SraRunInfo.${id}.txt sampleInfo.csv designTabFull.csv ${params.seqtype}
    echo \$(esearch --help 2>&1) | sed 's/^esearch //; s/Query .*\$//' > entrez_direct.version.txt
    """
}

