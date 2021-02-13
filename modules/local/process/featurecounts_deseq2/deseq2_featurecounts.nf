include { initOptions; saveFiles; getRealPath } from '../functions'

/*
 * STEP 7.4: Differential analysis with DESeq2
 */
process DESEQ2_FEATURECOUNTS {
    tag "${meta.antibody}"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/${meta.peaktype}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:"${meta.antibody}") }
        
    conda (params.conda ? "${params.conda_softwares.rbase}" : null)

    when:
    params.macs_gsize && !params.skip_consensus_peaks && !params.skip_diff_analysis

    input:
    tuple val(meta), path(counts)
    path deseq2_pca_header
    path deseq2_clustering_header
    val options

    output:
    path '*.tsv', emit: tsv
    path '*igv.txt', emit: igv
    path '*.{RData,results.txt,pdf,log}'
    path 'sizeFactors'
    path '*vs*/*.{pdf,txt}'
    path '*vs*/*.bed'

    script:
    prefix = "${meta.antibody}.consensus_peaks"
    vst = params.deseq2_vst ? '--vst TRUE' : ''
    def curr_path = getRealPath()
    def ioptions = initOptions(options)
    """
    ${curr_path}/utilities/install_packages.r optparse DESeq2 vsn ggplot2 RColorBrewer pheatmap lattice BiocParallel
    ${curr_path}/featurecounts_deseq2/featurecounts_deseq2.r \\
        --featurecount_file "$counts" \\
        --outdir ./ \\
        --outprefix "differential" \\
        --outsuffix '' \\
        --cores $task.cpus \\
        $vst

    sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
    sed -i -e 's/DESeq2 /${meta.antibody} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv

    sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
    sed -i -e 's/DESeq2 /${meta.antibody} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv

    find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "${meta.peaktype}/${ioptions.publish_dir}/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt
    """
}
