// Import generic module functions
include { saveFiles; getRealPath } from '../functions'

params.options = [:]

/*
 * Convert GTF file to BED format
 */
process GTF2BED {
    tag "$gtf"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', publish_id:'') }

    conda (params.conda ? "conda-forge::perl=5.26.2" : null)

    input:
    path gtf
    
    output:
    path '*.bed'

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def curr_path = getRealPath()
    """
    ${curr_path}/gtf2bed/gtf2bed $gtf > ${gtf.baseName}.bed
    """
}
