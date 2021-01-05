// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

/*
 * Consensus peaks across samples, create boolean filtering file, SAF file for featureCounts
 */
process MACS2_CONSENSUS {
    tag "$meta.id"
    label 'process_long'
    publishDir "${params.outdir}/${meta.peaktype}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.conda ? "${params.conda_softwares.bedtools} ${params.conda_softwares.rbase}" : null)

    input:
    tuple val(meta), path(peaks)
    val options

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.saf"), emit: saf
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.boolean.txt"), emit: boolean_txt
    tuple val(meta), path("*.intersect.txt"), emit: intersect_txt

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    if (meta.multiple_groups || meta.replicates_exist) {
        def software     = getSoftwareName(task.process)
        def ioptions     = initOptions(options)
        def prefix       = ioptions.suffix ? "${meta.id}${ioptions.suffix}.consensus_peaks" : "${meta.id}.consensus_peaks"
        def peak_type    = meta.peaktype
        def mergecols    = meta.peaktype=="narrowPeak" ? (2..10).join(',') : (2..9).join(',')
        def collapsecols = meta.peaktype=="narrowPeak" ? (['collapse']*9).join(',') : (['collapse']*8).join(',')
        def expandparam  = meta.peaktype=="narrowPeak" ? '--is_narrow_peak' : ''
        """
        sort -T '.' -k1,1 -k2,2n ${peaks.collect{it.toString()}.sort().join(' ')} \\
            | mergeBed -c $mergecols -o $collapsecols > ${prefix}.txt

        macs2_merged_expand.py \\
            ${prefix}.txt \\
            ${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks.${peak_type}","")} \\
            ${prefix}.boolean.txt \\
            --min_replicates $params.min_reps_consensus \\
            $expandparam

        awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed

        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
        awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${prefix}.boolean.txt >> ${prefix}.saf

        install_packages.r optparse UpSetR
        plot_peak_intersect.r -i ${prefix}.boolean.intersect.txt -o ${prefix}.boolean.intersect.plot.pdf
        """
    }
}
