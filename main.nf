#!/usr/bin/env nextflow
/*
========================================================================================
                         qiubio-nf-core/chipseq
========================================================================================
 qiubio-nf-core/chipseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/jianhong/chipseq
----------------------------------------------------------------------------------------
*/
if( !nextflow.version.matches('20.04+') ) {
    println "This workflow requires Nextflow version 20.04 or greater -- You are running version $nextflow.version"
    exit 1
}

nextflow.enable.dsl = 2

/*
 * Print help message if required
 */
if (params.help) {
    def command = "nextflow run jianhong/chipseq -r dev --input design.csv --genome GRCh37 -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help("$projectDir/nextflow_schema.json", command)
    exit 0
}

////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

/*
 * Reference genomes
 */
// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Configurable variables
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.gene_bed = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.macs_gsize = params.genome ? params.genomes[ params.genome ].macs_gsize ?: false : false
params.deep_gsize = params.genome ? params.genomes[ params.genome ].deep_gsize ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
anno_readme = params.genome ? params.genomes[ params.genome ].readme ?: false : false
params.species = params.genome ? params.genomes[ params.genome ].species ?: params.genome : params.species?:false

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

/*
 * Validate parameters
 */
if (params.input)     { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }
if (params.gtf)       { ch_gtf = file(params.gtf, checkIfExists: true) } else { exit 1, 'GTF annotation file not specified!' }
if (params.gene_bed)  { ch_gene_bed = file(params.gene_bed, checkIfExists: true) }
if (params.blacklist) { ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true) } else { ch_blacklist = Channel.empty() }

if (params.fasta) {
    ch_fasta = file(params.fasta, checkIfExists: true)
} else {
    exit 1, 'Fasta file not specified!'
}

// Save AWS IGenomes file containing annotation version
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// If --gtf is supplied along with --genome
// Make gene bed from supplied --gtf instead of using iGenomes one automatically
def makeBED = false
if (!params.gene_bed) {
    makeBED = true
} else if (params.genome && params.gtf) {
    if (params.genomes[ params.genome ].gtf != params.gtf) {
        makeBED = true
    }
}

/*
 * Check parameters
 */
Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles
Checks.macs2_warn(params, log)         // Show a big warning message if we're not running MACS

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

/*
 * Stage config files
 */
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

// JSON files required by BAMTools for alignment filtering
ch_bamtools_filter_se_config = file(params.bamtools_filter_se_config, checkIfExists: true)
ch_bamtools_filter_pe_config = file(params.bamtools_filter_pe_config, checkIfExists: true)

// Header files for MultiQC
ch_spp_nsc_header = file("$projectDir/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
ch_spp_rsc_header = file("$projectDir/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
ch_spp_correlation_header = file("$projectDir/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
ch_peak_count_header = file("$projectDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_annotation_header = file("$projectDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
ch_deseq2_pca_header = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_deseq2_clustering_header = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)

// deepTools genomic elements bed files
ch_genomic_elements_bed = params.genomicElements? Channel.fromPath(params.genomicElements, checkIfExists: true) : Channel.empty()

// index.Rmd
ch_index_docs = file("$projectDir/docs/index.Rmd", checkIfExists: true)

////////////////////////////////////////////////////
/* --          PARAMETER SUMMARY               -- */
////////////////////////////////////////////////////

/*
 * Print parameter summary
 */
// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
run_name = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    run_name = workflow.runName
}
summary = Schema.params_summary(workflow, params, run_name)
log.info Headers.nf_core(workflow, params.monochrome_logs)
log.info summary.collect { k,v -> "${k.padRight(22)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GTF2BED                             } from './modules/local/process/gtf2bed'
include { GET_CHROM_SIZES                     } from './modules/local/process/get_chrom_sizes'
include { MAKE_GENOME_FILTER                  } from './modules/local/process/make_genome_filter'
include { BEDTOOLS_GENOMECOV                  } from './modules/local/process/bedtools_genomecov'
include { PLOT_HOMER_ANNOTATEPEAKS            } from './modules/local/process/plot_homer_annotatepeaks'
include { PLOT_MACS2_QC                       } from './modules/local/process/plot_macs2_qc'
include { MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS } from './modules/local/process/multiqc_custom_phantompeakqualtools'
include { MULTIQC_CUSTOM_PEAKS                } from './modules/local/process/multiqc_custom_peaks'
include { MACS2_CONSENSUS                     } from './modules/local/process/macs2_consensus'
include { FRIP_SCORE                          } from './modules/local/process/frip_score'
//include { DESEQ2_FEATURECOUNTS                } from './modules/local/process/deseq2_featurecounts'
include { IGV                                 } from './modules/local/process/igv'
include { OUTPUT_DOCUMENTATION                } from './modules/local/process/doc/output_documentation'
include { GET_SOFTWARE_VERSIONS               } from './modules/local/process/get_software_versions'
include { MULTIQC                             } from './modules/local/process/multiqc'

include { INPUT_CHECK                         } from './modules/local/subworkflow/input_check'
include { BAM_CLEAN                           } from './modules/local/subworkflow/bam_clean'

include { JO_METAGENE_ANALYSIS                } from './modules/local/subworkflow/metagene_analysis'
include { JO_CHECKSUMS                        } from './modules/local/process/checksum/checksum'
include { JO_TRACKHUB                         } from './modules/local/process/ucsc_track/ucsc_track'
include { JO_INDEX                            } from './modules/local/process/create_index/create_index'
include { JO_DIFFBIND as JO_DIFFBIND_HOMER
          JO_DIFFBIND as JO_DIFFBIND_WITHOUT_CONTROL
          JO_DIFFBIND as JO_DIFFBIND_HOMER_WITHOUT_CONTROL
          JO_DIFFBIND                         } from './modules/local/process/diffbind/diffbind'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

include { BWA_INDEX                     } from './modules/nf-core/software/bwa/index/main'
include { PICARD_MERGESAMFILES          } from './modules/nf-core/software/picard/mergesamfiles/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from './modules/nf-core/software/picard/collectmultiplemetrics/main'
include { PRESEQ_LCEXTRAP               } from './modules/nf-core/software/preseq/lcextrap/main'
include { UCSC_BEDRAPHTOBIGWIG          } from './modules/nf-core/software/ucsc/bedgraphtobigwig/main'
include { DEEPTOOLS_COMPUTEMATRIX       } from './modules/nf-core/software/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE         } from './modules/nf-core/software/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP         } from './modules/nf-core/software/deeptools/plotheatmap/main'
include { DEEPTOOLS_PLOTFINGERPRINT     } from './modules/nf-core/software/deeptools/plotfingerprint/main'
include { PHANTOMPEAKQUALTOOLS          } from './modules/nf-core/software/phantompeakqualtools/main'
include { MACS2_CALLPEAK as MACS2_CALLPEAK_WITHOUT_CONTROL
          MACS2_CALLPEAK                } from './modules/nf-core/software/macs2/callpeak/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MACS2
          HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS } from './modules/nf-core/software/homer/annotatepeaks/main'
include { SUBREAD_FEATURECOUNTS         } from './modules/nf-core/software/subread/featurecounts/main'

include { FASTQC_TRIMGALORE             } from './modules/nf-core/subworkflow/fastqc_trimgalore'
include { MAP_BWA_MEM                   } from './modules/nf-core/subworkflow/map_bwa_mem'
include { MARK_DUPLICATES_PICARD        } from './modules/nf-core/subworkflow/mark_duplicates_picard'
include { HOMER_CALLPEAK as HOMER_CALLPEAK_WITHOUT_CONTROL
          HOMER_CALLPEAK                } from './modules/nf-core/subworkflow/findpeak_homer'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    /*
     * Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK (
        ch_input,
        params.seq_center,
        [:]
    )
    
    JO_CHECKSUMS (
        INPUT_CHECK.out.reads,
        [:]
    )

    /*
     * Prepare genome files
     */
    ch_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BWA_INDEX ( ch_fasta, params.modules['bwa_index'] ).index

    if (makeBED) { ch_gene_bed = GTF2BED ( ch_gtf, [:] ) }

    MAKE_GENOME_FILTER (
        GET_CHROM_SIZES ( ch_fasta, [:] ).sizes,
        ch_blacklist.ifEmpty([]),
        [:]
    )
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(MAKE_GENOME_FILTER.out.version.first().ifEmpty(null))

    /*
     * Read QC & trimming
     */
    nextseq = params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
    params.modules['trimgalore'].args += nextseq
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc,
        params.skip_trimming,
        params.modules['fastqc'],
        params.modules['trimgalore']
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    /*
     * Map reads & BAM QC
     */
    score = params.bwa_min_score ? " -T ${params.bwa_min_score}" : ''
    params.modules['bwa_mem'].args += score
    MAP_BWA_MEM (
        FASTQC_TRIMGALORE.out.reads,
        ch_index,
        ch_fasta,
        params.modules['bwa_mem'],
        params.modules['samtools_sort_lib']
    )
    ch_software_versions = ch_software_versions.mix(MAP_BWA_MEM.out.bwa_version.first())
    ch_software_versions = ch_software_versions.mix(MAP_BWA_MEM.out.samtools_version.first().ifEmpty(null))

    /*
     * Merge resequenced BAM files
     */
    MAP_BWA_MEM
        .out
        .bam
        .map {
            meta, bam ->
                fmeta = meta.findAll { !(it.key in meta.sample_uniq_keys) }
                fmeta.id = fmeta.id.split('_')[0..-2].join('_')
                [ fmeta, bam ] }
       .groupTuple(by: [0])
       .map { it ->  [ it[0], it[1].flatten() ] }
       .set { ch_sort_bam }

    PICARD_MERGESAMFILES (
        ch_sort_bam,
        params.modules['picard_mergesamfiles']
    )
    ch_software_versions = ch_software_versions.mix(PICARD_MERGESAMFILES.out.version.first().ifEmpty(null))

    /*
     * Mark duplicates & filter BAM files
     */
    MARK_DUPLICATES_PICARD (
        PICARD_MERGESAMFILES.out.bam,
        params.modules['picard_markduplicates'],
        params.modules['samtools_sort_merged_lib']
    )

    // Fix getting name sorted BAM here for PE/SE
    BAM_CLEAN (
        MARK_DUPLICATES_PICARD.out.bam.join(MARK_DUPLICATES_PICARD.out.bai, by: [0]),
        MAKE_GENOME_FILTER.out.bed.collect(),
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config,
        params.modules['bam_filter'],
        params.modules['bam_remove_orphans'],
        params.modules['samtools_sort_filter']
    )
    ch_software_versions = ch_software_versions.mix(BAM_CLEAN.out.bamtools_version.first().ifEmpty(null))

    /*
     * Post alignment QC
     */
    PICARD_COLLECTMULTIPLEMETRICS (
        BAM_CLEAN.out.bam,
        ch_fasta,
        params.modules['picard_collectmultiplemetrics']
    )

    PRESEQ_LCEXTRAP (
        BAM_CLEAN.out.bam,
        params.modules['preseq_lcextrap']
    )
    ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.version.first().ifEmpty(null))

    PHANTOMPEAKQUALTOOLS (
        BAM_CLEAN.out.bam,
        params.modules['phantompeakqualtools']
    )
    ch_software_versions = ch_software_versions.mix(PHANTOMPEAKQUALTOOLS.out.version.first().ifEmpty(null))

    MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS (
        PHANTOMPEAKQUALTOOLS.out.spp.join(PHANTOMPEAKQUALTOOLS.out.rdata, by: [0]),
        ch_spp_nsc_header,
        ch_spp_rsc_header,
        ch_spp_correlation_header,
        params.modules['multiqc_custom_phantompeakqualtools']
    )

    /*
     * Coverage tracks
     */
    BEDTOOLS_GENOMECOV (
        BAM_CLEAN.out.bam.join(BAM_CLEAN.out.flagstat, by: [0]),
        params.modules['bedtools_genomecov']
    )

    UCSC_BEDRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.bedgraph,
        GET_CHROM_SIZES.out.sizes,
        params.modules['ucsc_bedgraphtobigwig']
    )
    ch_software_versions = ch_software_versions.mix(UCSC_BEDRAPHTOBIGWIG.out.version.first().ifEmpty(null))

    /*
     * Coverage plots
     */
    DEEPTOOLS_COMPUTEMATRIX (
        UCSC_BEDRAPHTOBIGWIG.out.bigwig,
        ch_gene_bed,
        params.modules['deeptools_computematrix']
    )
    ch_software_versions = ch_software_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.version.first().ifEmpty(null))

    DEEPTOOLS_PLOTPROFILE (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix,
        params.modules['deeptools_plotprofile']
    )

    DEEPTOOLS_PLOTHEATMAP (
        DEEPTOOLS_COMPUTEMATRIX.out.matrix,
        params.modules['deeptools_plotheatmap']
    )

    /*
     * Refactor channels: [ val(meta), [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
     */
    BAM_CLEAN
        .out
        .bam
        .join ( BAM_CLEAN.out.bai, by: [0] )
        .map { meta, bam, bai -> meta.control ? null : [ meta.id, [ bam ] , [ bai ] ] }
        .set { ch_control_bam_bai }
        
    BAM_CLEAN
        .out
        .bam
        .join ( BAM_CLEAN.out.bai, by: [0] )
        .map { meta, bam, bai -> meta.control ? [ meta.control, meta, [ bam ], [ bai ] ] : null }
        .combine(ch_control_bam_bai, by: 0)
        .map { it -> [ it[1] , it[2] + it[4], it[3] + it[5] ] }
        .set { ch_ip_control_bam_bai }

    /*
     * plotFingerprint for IP and control together
     */
    params.modules['deeptools_plotfingerprint'].args += " --numberOfSamples $params.fingerprint_bins"
    DEEPTOOLS_PLOTFINGERPRINT (
        ch_ip_control_bam_bai,
        params.modules['deeptools_plotfingerprint']
    )

    
    /*
     * do metagene analysis for given bed files
     */
    JO_METAGENE_ANALYSIS(
        BAM_CLEAN.out.bam,
        UCSC_BEDRAPHTOBIGWIG.out.bigwig,
        ch_genomic_elements_bed,
        params.modules['jo_merge_rep_bam'],
        params.modules['jo_metagene']
    )
    
    if (params.macs_gsize) {

        /*
         * Call peaks
         */
        pileup = params.save_macs_pileup ? '--bdg --SPMR' : ''
        fdr = params.macs_fdr ? "--qvalue ${params.macs_fdr}" : ''
        pvalue = params.macs_pvalue ? "--pvalue ${params.macs_pvalue}" : ''
        params.modules['macs2_callpeak'].args += " $pileup $fdr $pvalue"

        // call peaks without input
        BAM_CLEAN.out.bam.map{ meta, bam -> [meta, bam, []]}
                         .set{ch_ip_bam_no_ctl}

        MACS2_CALLPEAK_WITHOUT_CONTROL (
            ch_ip_bam_no_ctl,
            params.macs_gsize,
            params.modules['macs2_callpeak_without_control']
        )
        
        BAM_CLEAN.out.bam
            .join(MACS2_CALLPEAK_WITHOUT_CONTROL.out.peak, by: [0])
            .set { ch_ip_peak_no_control }
        ch_ip_peak_no_control
            .map{meta, bam, peak -> [meta.antibody, meta]}
            .groupTuple()
            .join(ch_ip_peak_no_control.map{[it[0].antibody, it[1]]}.groupTuple()
                    .join(ch_ip_peak_no_control.map{[it[0].antibody, it[2]]}.groupTuple()))
            .map{[it[1], it[2], it[3]]}
            .set{ch_diffbind_without_control}
        JO_DIFFBIND_WITHOUT_CONTROL (
            ch_diffbind_without_control,
            ch_gtf,
            ch_blacklist.ifEmpty([]),
            params.modules['jo_diffbind_macs2_without_control']
        )
        
        // Create channel: [ val(meta), ip_bam, control_bam ]
        ch_ip_control_bam_bai
            .map { meta, bams, bais -> [ meta , bams[0], bams[1] ] }
            .set { ch_ip_control_bam }

        MACS2_CALLPEAK (
            ch_ip_control_bam,
            params.macs_gsize,
            params.modules['macs2_callpeak']
        )
        ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK.out.version.first().ifEmpty(null))

        ch_ip_control_bam
            .join(MACS2_CALLPEAK.out.peak, by: [0])
            .map { it -> [ it[0], it[1], it[3] ] }
            .set { ch_ip_peak }
        FRIP_SCORE (
            ch_ip_peak,
            params.modules['frip_score']
        )

        ch_ip_peak
            .join(FRIP_SCORE.out.txt, by: [0])
            .map { it -> [ it[0], it[2], it[3] ] }
            .set { ch_ip_peak_frip }
        MULTIQC_CUSTOM_PEAKS (
            ch_ip_peak_frip,
            ch_peak_count_header,
            ch_frip_score_header,
            params.modules['multiqc_custom_peaks']
        )

        params.modules['plot_macs2_qc'].publish_dir += "/macs2_qc"
        PLOT_MACS2_QC (
            MACS2_CALLPEAK.out.peak.map{[it[0].peaktype, it[1]]}.groupTuple(),
            params.modules['plot_macs2_qc']
        )

        HOMER_ANNOTATEPEAKS_MACS2 (
            MACS2_CALLPEAK.out.peak,
            ch_fasta,
            ch_gtf,
            params.modules['homer_annotatepeaks_macs2'],
            false
        )
        ch_software_versions = ch_software_versions.mix(HOMER_ANNOTATEPEAKS_MACS2.out.version.first().ifEmpty(null))

        params.modules['plot_homer_annotatepeaks'].publish_dir += "/homer_anno_qc"
        PLOT_HOMER_ANNOTATEPEAKS (
            HOMER_ANNOTATEPEAKS_MACS2.out.txt.map{[it[0].peaktype, it[1]]}.groupTuple(),
            ch_peak_annotation_header,
            "_peaks.annotatePeaks.txt",
            params.modules['plot_homer_annotatepeaks']
        )

        // Create channel: [ meta , [ peaks ] ]
        // Where meta = [ id:antibody, multiple_groups:true/false, replicates_exist:true/false ]
        MACS2_CALLPEAK
            .out
            .peak
            .map { meta, peak -> [ meta.antibody, meta.id.split('_')[0..-2].join('_'), meta.peaktype, peak ] }
            .groupTuple()
            .map {
                antibody, groups, peaktype, peaks ->
                    [ antibody,
                      groups.groupBy().collectEntries { [(it.key) : it.value.size()] },
                      peaktype[0],
                      peaks ] }
            .map {
                antibody, groups, peaktype, peaks ->
                    def meta = [:]
                    meta.id = antibody
                    meta.peaktype = peaktype
                    meta.multiple_groups = groups.size() > 1
                    meta.replicates_exist = groups.max { groups.value }.value > 1
                    [ meta, peaks ] }
            .set { ch_antibody_peaks }

        params.modules['macs2_consensus'].publish_dir += "/consensus"
        MACS2_CONSENSUS (
            ch_antibody_peaks,
            params.modules['macs2_consensus']
        )

        params.modules['homer_annotatepeaks_consensus'].publish_dir += "/consensus"
        HOMER_ANNOTATEPEAKS_CONSENSUS (
            MACS2_CONSENSUS.out.bed,
            ch_fasta,
            ch_gtf,
            params.modules['homer_annotatepeaks_consensus'],
            false
        )
        // cut -f2- ${prefix}.annotatePeaks.txt | awk 'NR==1; NR > 1 {print \$0 | "sort -T '.' -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
        // paste $bool tmp.txt > ${prefix}.boolean.annotatePeaks.txt

        // Create channel: [ val(meta), ip_bam ]
        MACS2_CONSENSUS
            .out
            .saf
            .map { meta, saf -> [ meta.id, meta, saf ] }
            .set { ch_ip_saf }
        //ch_ip_saf.view()

        ch_ip_control_bam
            .map { meta, ip_bam, control_bam -> [ meta.antibody, meta, ip_bam ] }
            .combine(ch_ip_saf)
            .filter {
                antibody, meta, bam, ab, meta2, saf ->
                   ab == meta.antibody
            }
            .map {
                it ->
                    fmeta = it[1]
                    fmeta['replicates_exist'] = it[4]['replicates_exist']
                    fmeta['multiple_groups'] = it[4]['multiple_groups']
                    [ fmeta, it[2], it[5] ] }
            .set { ch_ip_bam }
        
        //ch_ip_bam.view()

        params.modules['subread_featurecounts'].publish_dir += "/consensus"
        SUBREAD_FEATURECOUNTS (
            ch_ip_bam,
            params.modules['subread_featurecounts']
        )
        ch_software_versions = ch_software_versions.mix(SUBREAD_FEATURECOUNTS.out.version.first().ifEmpty(null))

        //
        // DESEQ2_FEATURECOUNTS (
        //     params.modules['deseq2_featurecounts']
        // )
        // ch_deseq2_pca_header = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
        // ch_deseq2_clustering_header = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)

        ch_ip_peak
            .map{meta, bam, peak -> [meta.antibody, meta]}
            .groupTuple()
            .join(ch_ip_peak.map{[it[0].antibody, it[1]]}.groupTuple()
                    .join(ch_ip_peak.map{[it[0].antibody, it[2]]}.groupTuple()))
            .map{[it[1], it[2], it[3]]}
            .set{ch_diffbind}
        
        JO_DIFFBIND (
            ch_diffbind,
            ch_gtf,
            ch_blacklist.ifEmpty([]),
            params.modules['jo_diffbind_macs2']
        )
    }

    if(params.homer){
        /*
         * Call peaks
         */
        // call peaks without input
        BAM_CLEAN.out.bam.map{ meta, bam -> [meta, bam, []]}
                         .set{ch_ip_bam_no_ctl_homer}

        HOMER_CALLPEAK_WITHOUT_CONTROL (
            ch_ip_bam_no_ctl_homer,
            ch_fasta,
            ch_gtf,
            params.modules['homer_maketagdirecotry'],
            params.modules['homer_findpeaks'],
            params.modules['homer_annotatepeaks'],
            params.modules['homer_pos2bed'],
            [publish_dir:"peaks_without_control"]
        )
        BAM_CLEAN.out.bam
            .join(HOMER_CALLPEAK_WITHOUT_CONTROL.out.bed, by: [0])
            .set { ch_ip_peak_homer_no_control }
        ch_ip_peak_homer_no_control
            .map{meta, bam, peak -> [meta.antibody, meta]}
            .groupTuple()
            .join(ch_ip_peak_homer_no_control.map{[it[0].antibody, it[1]]}.groupTuple()
                    .join(ch_ip_peak_homer_no_control.map{[it[0].antibody, it[2]]}.groupTuple()))
            .map{[it[1], it[2], it[3]]}
            .set{ch_diffbind_homer_without_control}
        JO_DIFFBIND_HOMER_WITHOUT_CONTROL (
            ch_diffbind_homer_without_control,
            ch_gtf,
            ch_blacklist.ifEmpty([]),
            params.modules['jo_diffbind_homer_without_control']
        )
        
        // Create channel: [ val(meta), ip_bam, control_bam ]
        ch_ip_control_bam_bai
            .map { meta, bams, bais -> [ meta , bams[0], bams[1] ] }
            .set { ch_ip_control_bam_homer }
            
        HOMER_CALLPEAK (
            ch_ip_control_bam_homer,
            ch_fasta,
            ch_gtf,
            params.modules['homer_maketagdirecotry'],
            params.modules['homer_findpeaks'],
            params.modules['homer_annotatepeaks'],
            params.modules['homer_pos2bed'],
            [publish_dir:"peaks_with_control"]
        )
        
        ch_ip_control_bam_homer
            .map { it -> [it[0].id, it[0], it[1]]}
            .join(HOMER_CALLPEAK.out.bed.map{[it[0].id, it[0], it[1]]}, by: [0])
            .map { it -> [ it[1], it[2], it[4] ] }
            .set { ch_ip_peak_homer }

        ch_ip_peak_homer
            .map{meta, bam, peak -> [meta.antibody, meta]}
            .groupTuple()
            .join(ch_ip_peak_homer.map{[it[0].antibody, it[1]]}.groupTuple()
                    .join(ch_ip_peak_homer.map{[it[0].antibody, it[2]]}.groupTuple()))
            .map{[it[1], it[2], it[3]]}
            .set{ch_diffbind_homer}
        
        JO_DIFFBIND_HOMER (
            ch_diffbind_homer,
            ch_gtf,
            ch_blacklist.ifEmpty([]),
            params.modules['jo_diffbind_homer']
        )
    
    }
    /*
     * Create IGV session
     */
    IGV (
        ch_fasta,
        UCSC_BEDRAPHTOBIGWIG.out.bigwig.collect{it[1]}.ifEmpty([]),
        MACS2_CALLPEAK.out.peak.collect{it[1]}.ifEmpty([]),
        MACS2_CONSENSUS.out.bed.collect{it[1]}.ifEmpty([]),
        params.modules['ucsc_bedgraphtobigwig'],
        params.modules['macs2_callpeak'],
        params.modules['macs2_consensus'],
        [:]
    )
    
    /*
     * Create ucsc trackhub
     */
   
   JO_METAGENE_ANALYSIS.out.bw.map{[it[0], it[1]]}.collect().ifEmpty([])
        .concat(UCSC_BEDRAPHTOBIGWIG.out.bigwig.map{[it[0], it[1]]}.collect().ifEmpty([]), 
                MACS2_CALLPEAK.out.peak.map{[it[0], it[1]]}.collect().ifEmpty([]),
                MACS2_CONSENSUS.out.bed.map{[it[0], it[1]]}.collect().ifEmpty([]),
                HOMER_CALLPEAK.out.bed.map{[it[0], it[1]]}.collect().ifEmpty([]))
        .flatten()
        .collate(2)
        .set{ch_trackhub}
   ch_trackhub.map{[it[0].id]}.collect().set{ch_trackhub_name}
   ch_trackhub.map{[it[1]]}.collect().set{ch_trackhub_track}

   JO_TRACKHUB(
        ch_trackhub_name,
        ch_trackhub_track,
        ch_input,
        GET_CHROM_SIZES.out.sizes,
        params.modules['jo_trackhub']
   )

    /*
     * Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect(),
        params.modules['get_software_versions']
    )

    OUTPUT_DOCUMENTATION (
        ch_output_docs,
        ch_output_docs_images,
        [:]
    )

    /*
     * MultiQC
     */
    workflow_summary = Schema.params_mqc_summary(summary)
    ch_workflow_summary = Channel.value(workflow_summary)
    params.modules['multiqc'].publish_dir += "/$run_name"
    MULTIQC (
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        GET_SOFTWARE_VERSIONS.out.yaml.collect(),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),

        FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
        FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
        FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),

        MAP_BWA_MEM.out.stats.collect{it[1]},
        MAP_BWA_MEM.out.flagstat.collect{it[1]},
        MAP_BWA_MEM.out.idxstats.collect{it[1]},

        MARK_DUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]),
        MARK_DUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]),
        MARK_DUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]),
        MARK_DUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]),

        BAM_CLEAN.out.stats.collect{it[1]}.ifEmpty([]),
        BAM_CLEAN.out.flagstat.collect{it[1]}.ifEmpty([]),
        BAM_CLEAN.out.idxstats.collect{it[1]}.ifEmpty([]),
        PICARD_COLLECTMULTIPLEMETRICS.out.metrics.collect{it[1]}.ifEmpty([]),

        PRESEQ_LCEXTRAP.out.ccurve.collect{it[1]}.ifEmpty([]),
        DEEPTOOLS_PLOTPROFILE.out.table.collect{it[1]}.ifEmpty([]),
        DEEPTOOLS_PLOTFINGERPRINT.out.matrix.collect{it[1]}.ifEmpty([]),
        PHANTOMPEAKQUALTOOLS.out.spp.collect{it[1]}.ifEmpty([]),
        MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.nsc.collect{it[1]}.ifEmpty([]),
        MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.rsc.collect{it[1]}.ifEmpty([]),
        MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.correlation.collect{it[1]}.ifEmpty([]),

        MULTIQC_CUSTOM_PEAKS.out.count.collect{it[1]}.ifEmpty([]),
        MULTIQC_CUSTOM_PEAKS.out.frip.collect{it[1]}.ifEmpty([]),
        PLOT_HOMER_ANNOTATEPEAKS.out.tsv.collect().ifEmpty([]),
        SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]),
        // path ('macs/consensus/*') from ch_macs_consensus_deseq_mqc.collect().ifEmpty([])

        params.modules['multiqc']
    )
    
    
    JO_INDEX (
        ch_index_docs,
        ch_input,
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
        MULTIQC.out.plots,
        ch_output_docs_images,
        JO_CHECKSUMS.out.md5.collect(),
        GET_SOFTWARE_VERSIONS.out.yaml.collect(),
        [:]
    )
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    def multiqc_report = []
    Completion.email(workflow, params, summary, run_name, projectDir, multiqc_report, log)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////