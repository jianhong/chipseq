// Import generic module functions
include { saveFiles } from './functions'

/*
 * Create IGV session file
 */
process IGV {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:'') }

    input:
    path fasta
    path ("${bigwig_options.publish_dir}/*")
    path ("${peak_options.publish_dir}/*")
    path ("${consensus_options.publish_dir}/*")
    val bigwig_options
    val peak_options
    val consensus_options
    val options

    output:
    path "*files.txt", emit: txt
    path "*.xml", emit: xml

    script: // scripts are bundled with the pipeline in nf-core/chipseq/bin/
    """
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; > bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; > peaks.igv.txt

    cat *.txt > igv_files.txt
    igv_files_to_session.py igv_session.xml igv_files.txt ../../genome/${fasta.getName()} --path_prefix '../../'
    """
}
