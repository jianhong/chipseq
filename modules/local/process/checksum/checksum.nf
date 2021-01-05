// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * MD5 checksums
 */
process JO_CHECKSUMS {
    tag "$name"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.tokenize(':')[-1].tokenize('_')[1].toLowerCase(), publish_id:name) }
    
    input:
    tuple val(meta), path(reads)
    val options

    output:
    path "md5.*.txt", emit: md5
    
    script:
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}.${ioptions.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        touch md5.${prefix}.txt
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        gunzip -c ${prefix}.fastq.gz > ${prefix}.fastq
        ${params.md5sum} ${prefix}.fastq >>md5.${prefix}.txt
        if [ "${meta.md5_1}" != "null" ]; then
            md5=(\$(${params.md5sum} ${prefix}.fastq.gz))
            if [ "\$md5" != "${meta.md5_1}" ]
            then
                echo "${meta.id} has checksum ${meta.md5_1}, but we got checksum \$md5!"
                exit 128
            fi
        fi
        """
    } else {
        """
        touch md5.${prefix}.txt
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        gunzip -c ${prefix}_1.fastq.gz > ${prefix}_1.fastq
        ${params.md5sum} ${prefix}_1.fastq >>md5.${prefix}.txt
        gunzip -c ${prefix}_2.fastq.gz > ${prefix}_2.fastq
        ${params.md5sum} ${prefix}_2.fastq >>md5.${prefix}.txt
        if [ "${meta.md5_1}" != "null" ]; then
            md5=(\$(${params.md5sum} ${prefix}_1.fastq.gz))
            if [ "\$md5" != "${meta.md5_1}" ]
            then
                echo "${meta.id} has checksum ${meta.md5_1}, but we got checksum \$md5!"
                exit 128
            fi
        fi
        if [ "${meta.md5_2}" != "null" ]; then
            md5=(\$(${params.md5sum} ${prefix}_2.fastq.gz))
            if [ "\$md5" != "${meta.md5_2}" ]
            then
                echo "${meta.id} has checksum ${meta.md5_2}, but we got checksum \$md5!"
                exit 128
            fi
        fi
        """
    }
}