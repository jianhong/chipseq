// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

def VERSION = '4.11'

process HOMER_ANNOTATEPEAKS {
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
    path fasta
    path gtf
    val options

    output:
    tuple val(meta), path("*annotatePeaks.txt"), emit: txt
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    annotatePeaks.pl \\
        $peak \\
        $fasta \\
        $ioptions.args \\
        -gtf $gtf \\
        -cpu $task.cpus \\
        > ${prefix}.annotatePeaks.txt

    echo $VERSION > ${software}.version.txt
    """
}
