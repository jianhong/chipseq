// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * Create trackhub
 */
process JO_TRACKHUB {
    publishDir "${params.outdir}",
        mode: 'copyNoFollow',
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:'', publish_id:'') }

    conda (params.conda ? "${params.modules_dir}/ucsc_track/environment.txt" : null)

    when:
    !params.skip_trackhub
    
    input:
    val name
    path bw
    path designtab
    val options

    output:
    path "trackhub/*"
    path "*.version.txt", emit: version
    
    script:
    def sampleLabel = name.join('___')
    def bws         = bw.join('___')
    """
    find * -type l -name "*.bigWig" -exec echo {} \\; > bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo {} \\; > peaks.igv.txt

    cat *.txt > track_files.txt
    
    create_trackhub.py trackhub track_files.txt $params.species $params.email $designtab '.mLb.clN__.mRp.clN' --path_prefix '../../../../'
    
    python -c "import trackhub; print(trackhub.__version__)" > trackhub.version.txt
    """
}