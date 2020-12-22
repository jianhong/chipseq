// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * Metagene analysis
 */
process JO_METAGENE {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.toLowerCase(), publish_id:name) }

    container "quay.io/biocontainers/deeptools:3.4.3--py_0"
    conda (params.conda ? "${params.conda_softwares.deeptools}" : null)

    when:
    params.genomicElements
    
    input:
    tuple val(name), path(bw)
    path bed
    val options

    output:
    tuple val(name), path("*.gz"), emit: count
    tuple val(name), path("*.pdf"), emit: pdf
    tuple val(name), path("*.tab"), emit: mat
    path "*.version.txt", emit: version
    
    script:
    def bigwig       = bw.join(' ')
    sampleLabel      = name.join(' ')
    """
    computeMatrix scale-regions \\
        --regionsFileName $bed \\
        --scoreFileName ${bigwig} \\
        --samplesLabel ${sampleLabel} \\
        --outFileName ${bed.getSimpleName()}.scale_regions.mat.gz \\
        --outFileNameMatrix ${bed.getSimpleName()}.scale_regions.vals.mat.tab \\
        --regionBodyLength $params.deepToolsBodySize \\
        --beforeRegionStartLength $params.deepToolsRegionSize \\
        --afterRegionStartLength $params.deepToolsRegionSize \\
        --skipZeros \\
        --numberOfProcessors $task.cpus
    plotProfile --matrixFile ${bed.getSimpleName()}.scale_regions.mat.gz \\
        --outFileName ${bed.getSimpleName()}.scale_regionsProfile.pdf \\
        --outFileNameData ${bed.getSimpleName()}.scale_regionsProfile.tab
    plotHeatmap --matrixFile ${bed.getSimpleName()}.scale_regions.mat.gz \\
        --outFileName ${bed.getSimpleName()}.scale_regionsHeatmap.pdf \\
        --outFileNameMatrix ${bed.getSimpleName()}.scale_regionsHeatmap.mat.tab
        
    computeMatrix reference-point \\
        --regionsFileName $bed \\
        --scoreFileName ${bigwig} \\
        --samplesLabel ${sampleLabel} \\
        --outFileName ${bed.getSimpleName()}.reference_${params.deepToolsReferencePoint}.mat.gz \\
        --outFileNameMatrix ${bed.getSimpleName()}.reference_${params.deepToolsReferencePoint}.vals.mat.tab \\
        --referencePoint $params.deepToolsReferencePoint \\
        --beforeRegionStartLength $params.deepToolsRegionSize \\
        --afterRegionStartLength $params.deepToolsRegionSize \\
        --skipZeros \\
        --numberOfProcessors $task.cpus
    plotProfile --matrixFile ${bed.getSimpleName()}.reference_${params.deepToolsReferencePoint}.mat.gz \\
        --outFileName ${bed.getSimpleName()}.reference_${params.deepToolsReferencePoint}Profile.pdf \\
        --outFileNameData ${bed.getSimpleName()}.reference_${params.deepToolsReferencePoint}Profile.tab
    plotHeatmap --matrixFile ${bed.getSimpleName()}.reference_${params.deepToolsReferencePoint}.mat.gz \\
        --outFileName ${bed.getSimpleName()}.reference_${params.deepToolsReferencePoint}Heatmap.pdf \\
        --outFileNameMatrix ${bed.getSimpleName()}.reference_${params.deepToolsReferencePoint}Heatmap.mat.tab
    
    computeMatrix --version | sed -e "s/computeMatrix //g" > deeptools.version.txt
    """
}