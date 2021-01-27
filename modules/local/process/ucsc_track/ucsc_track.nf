// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * Create trackhub
 */
process JO_TRACKHUB {
    label 'error_ignore'
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
    path size
    val options

    output:
    path "trackhub/*"
    path "*.version.txt", emit: version
    
    script:
    def sampleLabel = name.join(' ')
    def bws         = bw.join(' ')
    """
    n=($sampleLabel)
    s=($bws)
    paste <(printf '%s\n' "\${n[@]}") <(printf '%s\n' "\${s[@]}") > track_files.txt
    
    create_trackhub.py track_files.txt $params.species $size $params.email $designtab
    
    rsync trackhub/${params.species} tmp -a --copy-links -v
    rm -r trackhub/${params.species}
    mv tmp/${params.species} trackhub/
    
    python -c "import trackhub; print(trackhub.__version__)" > trackhub.version.txt
    """
}