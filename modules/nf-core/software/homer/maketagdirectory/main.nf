// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

def VERSION = '4.11'

process HOMER_MAKETAGDIRECTORY {
    tag "$meta.id"
    label 'process_medium'

    container (params.universalContainer? "${params.container}":"quay.io/biocontainers/homer:4.11--pl526h9a982cc_2")
    //container "https://depot.galaxyproject.org/singularity/homer:4.11--pl526h9a982cc_2"

    conda (params.conda ? "${params.conda_softwares.homer}" : null)

    input:
    tuple val(meta), path(bam), path(bamcontrol)
    val options

    output:
    tuple val(meta), path("${meta.id}_Tagdir"), path("${meta.id}_controlTagdir"), emit: tagdir
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def pe       = meta.single_end ? "":"-sspe"
    if (bamcontrol){
        """
        makeTagDirectory ${meta.id}_Tagdir $pe ${bam}
        makeTagDirectory ${meta.id}_controlTagdir $pe ${bamcontrol}
    
        echo $VERSION > ${software}.version.txt
        """
    } else {
        """
        makeTagDirectory ${meta.id}_Tagdir $pe ${bam}
        touch ${meta.id}_controlTagdir
    
        echo $VERSION > ${software}.version.txt
        """
    }
}
