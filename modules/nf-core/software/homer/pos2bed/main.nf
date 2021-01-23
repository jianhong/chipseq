// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

def VERSION = '4.11'

process HOMER_POS2BED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.peaktype}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container (params.universalContainer? "${params.container}":"quay.io/biocontainers/homer:4.11--pl526h9a982cc_2")
    //container "https://depot.galaxyproject.org/singularity/homer:4.11--pl526h9a982cc_2"

    conda (params.conda ? "${params.conda_softwares.homer}" : null)

    input:
    tuple val(meta), path(peak)
    val options

    output:
    tuple val(meta), path("${meta.id}_homer_${meta.peaktype}.bed"), emit: bed
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    pos2bed ${peak} -o ${meta.id}_homer_${meta.peaktype}.bed -track ${meta.id}

    echo $VERSION > ${software}.version.txt
    """
}
