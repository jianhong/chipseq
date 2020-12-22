// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

def VERSION = '1.2.2'

process PHANTOMPEAKQUALTOOLS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/phantompeakqualtools:1.2.2--0"
    //container "https://depot.galaxyproject.org/singularity/phantompeakqualtools:1.2.2--0"

    conda (params.conda ? "${params.conda_softwares.rbase}" : null)

    input:
    tuple val(meta), path(bam)
    val options

    output:
    tuple val(meta), path("*.out"), emit: spp
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.Rdata"), emit: rdata
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    install_packages.r snow snowfall bitops caTools Rsamtools hms-dbmi/spp
    run_spp.R -c="$bam" -savp="${prefix}.spp.pdf" -savd="${prefix}.spp.Rdata" -out="${prefix}.spp.out" -p=$task.cpus
    echo $VERSION > ${software}.version.txt
    """
}
