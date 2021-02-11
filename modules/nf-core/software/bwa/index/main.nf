// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process BWA_INDEX {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container (params.universalContainer? "${params.container}":"biocontainers/bwa:v0.7.17_cv1")
    //container "https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7"

    conda (params.conda ? "${params.conda_softwares.bwa}" : null)

    input:
    path fasta

    output:
    path "${fasta}.*", emit: index
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(params.options)
    """
    bwa index $ioptions.args $fasta
    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
    """
}
